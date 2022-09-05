//
//  FES.cpp
//  LosTopos
//
//  Created by Sadashige ISHIDA on 2018/05/06.
//
#ifdef FEOS

#include <stdio.h>
#include <numeric>
#include <exception>
#include <limits>
#include <unordered_map>
#include <array>
#include <math.h>

#include <Eigen/Geometry>
#include <Eigen/SparseCore>
#include <Eigen/SparseCholesky>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseLU>

#include <igl/cotmatrix.h>
#include <igl/invert_diag.h>
#include <igl/volume.h>
#include <igl/massmatrix.h>
#include <igl/per_vertex_normals.h>
#include <igl/principal_curvature.h>

//#include <igl/boundary_loop.h>

#include "Options.h"
#include "HGF.h"
#include "nondestructivetrimesh.h"

//#define advect_experiment
//#pragma optimize("", off)

double HGF::project_tri(double dt, Equation equation){
    switch (equation) {
        case EulerEquation:
            project_tri(dt);
            break;
        case SoapFilmEquations:
            project_tri(dt, uthvv(), thvv());
            break;
    }
    return dt;
}

double HGF::project_tri( double dt, const Eigen::VectorXd& Uth, const Eigen::VectorXd& Th){
    
    computePressureOnVertices(dt,Uth,Th);
    
    vert2triGrad( P_vert, GradP_tri );
    
    const size_t num_t=mesh().nt();
    
    if(Th.size()>0){
        Eigen::VectorXd Th_tri;
        scalarVert2tri(Th, Th_tri);
        for( size_t tti = 0; tti < real_triangles.size(); ++tti ){
            GradP_tri.row( tti )/=Th_tri[tti];
        }
    }
    
    for( size_t tti = 0; tti < real_triangles.size(); ++tti ){
        size_t ti=real_triangles[tti];
        
        triangle_velocities[ ti ] -= GradP_tri.row( ti );
        projectToTangent( triangle_velocities[ ti ], triangle_normals[ ti ], false );
    }
    
    return dt;
}

double HGF::step_SoapFilmEquations( double dt ){
    
    using namespace Eigen;
    size_t num_v = mesh().nv();
    size_t num_e=mesh().ne();
    size_t num_t=mesh().nt();
    
    if(with_gravity){
        addGravity2Flow(dt);
    }
    if(with_evaporation){
        evaporate(dt);
    }
    
    project_flow_v_tangent();
    
    set_flow_boundary();
    
    const std::vector<Vector3d> old_fv=fvvec();
    
    // Migrate Th (thickness) from vertices to a vector.
    // Migrate Uth (the time derivative of thickness) from vertices to a vector .
    
    VectorXd Th = thvv();
    VectorXd Uth = uthvv();
    const VectorXd PrevTh=Th;
    
    // DUth is flow_tension * H n.
    // We are using two surface tension related coefficients.
    // "surface_tension" is for vertical acceleration.
    // "flow_tension" is for thickness update.
    VectorXd DUth;
    computeSurfaceTensionForce(Th,DUth);
    
    // Compute Uth+=sigma DUth
    const double thickness_band=max_th_nm-min_th_nm;
    
    Uth+=dt*flow_tension*DUth;
    Th+=dt*Uth;
    clamp(Th,min_th_nm,max_th_nm);
    Uth=(Th-PrevTh)/dt;
    
    set_flow_boundary();
    
    // Project
    project_tri(dt, Uth, Th);
    
    // Advect velocity
    advect_triangleVelocity(dt,triangleAdvectType);
    
    // Advect Th and Uth
    advect_vertex(dt, Th,  vertexAdvectType);
    advect_vertex(dt, Uth,  vertexAdvectType);
    
    if(vorticityConfinementCoef>0){
        MatrixXd TriVelocity=getTriVelocity();
        vorticityConfinement(TriVelocity, vorticityConfinementCoef);
        setTriVelocity(TriVelocity);
    }
    
    correctLiquidVolumes(dt, Th, Uth);
    
    // Restore thickness and the time derivative to vertices.
    setThv(Th);
    setUThv(Uth);
    
    return dt;
}

void HGF::set_flow_boundary(){
    for( size_t gti = 0; gti<ghost_triangles.size();++gti){
        triangle_velocities[ghost_triangles[gti]].setZero();
    }
    
    for( size_t gvi = 0; gvi<ghost_vertices.size();++gvi){
        thv(ghost_vertices[gvi])=default_th;
        uthv(ghost_vertices[gvi])=0.0;
        fv(ghost_vertices[gvi]).setZero();
    }
    
    // Set thv and uth for boundary vertices
    const size_t num_v=mesh().nv();
    
    for(size_t vi=0;vi<num_v;++vi){
        // Apply only boundary vertices
        // Ignore non_solid vertices and ghost vertices
        if( not surfTrack()->vertex_is_any_solid(vi) or verts_to_real[vi]==-1 ){
            continue;
        }
        
        // Take the nearest value from non-ghost and non-boundary neighbors
        
        std::vector<size_t> adjacent_vertices;
        mesh().get_adjacent_vertices(vi, adjacent_vertices);
        const size_t num_adj=adjacent_vertices.size();
        
        std::vector<double> weights(num_adj);
        std::vector<double> distances(num_adj);
        
        double minSquareDist=1e8;
        int nearestVertex=-1;
        for(size_t adj_v:adjacent_vertices){
            if( surfTrack()->vertex_is_any_solid( adj_v ) ){
                continue;
            }
            const double squareDistance=(pos(vi)-pos(adj_v)).norm();
            if(squareDistance<minSquareDist) minSquareDist=squareDistance;
            nearestVertex=adj_v;
        }
        assert(nearestVertex!=-1);
        
        thv(vi)=thv(nearestVertex);
        uthv(vi)=uthv(nearestVertex);
        
    }
    
    for(size_t rti=0;rti<real_triangles.size();++rti){
        const size_t ti = real_triangles[rti];
        const auto& edge_inds=mesh().m_triangle_to_edge_map[ti];
        
        for(size_t ei_t=0;ei_t<3;++ei_t ){
            
            const size_t edge_ind=edge_inds[ei_t];
            if(surfTrack()->edge_is_all_solid(edge_ind)){
                const auto&edge_verts=mesh().m_edges[edge_ind];
                
                const Eigen::Vector3d edge=pos(edge_verts[1])-pos(edge_verts[0]);
                
            }
        }
    }
    for(size_t gti=0;gti<ghost_triangles.size();++gti){
        const size_t ti = ghost_triangles[gti];
        triangle_velocities[ti]<<0.0,0.0,0.0;
    }
    
}

double HGF::step_EulerEquation( double dt ){
    
    // set number of simplices
    using namespace Eigen;
    size_t num_v = mesh().nv();
    size_t num_e=mesh().ne();
    size_t num_t=mesh().nt();
    
    if(with_gravity){
        addGravity2Flow(dt);
    }
    
    if(with_evaporation){
        evaporate(dt);
    }
    
    project_flow_v_tangent();
    
    set_flow_boundary();
    
    project_tri(dt);
    
    // advect velocity
    advect_triangleVelocity( dt,triangleAdvectType);
    
    VectorXd Th=thvv();
    
    // advect scalar values
    advect_vertex( dt, Th, vertexAdvectType );
    clamp(Th,min_th_nm,max_th_nm);
    setThv(Th);
    
    if(vorticityConfinementCoef>0){
        MatrixXd TriVelocity=getTriVelocity();
        vorticityConfinement(TriVelocity, vorticityConfinementCoef);
        setTriVelocity(TriVelocity);
    }
    
    project_tri(dt);
    
    set_flow_boundary();
    
    return actual_dt;
}

// Thif function attempts to use triangle thickness and mass preserving advection on it.
double HGF::step_EulerEquation_TriangleTh(double dt){
    // set number of simplices
    using namespace Eigen;
    size_t num_v = mesh().nv();
    size_t num_e=mesh().ne();
    size_t num_t=mesh().nt();
    
    static bool once=false;
    if(not once){
        VectorXd Th_tri_init;
        scalarVert2tri(thvv(), Th_tri_init);
        triangle_thickness=std::vector<double>(Th_tri_init.data(),Th_tri_init.data()+Th_tri_init.size());
        once=true;
    }
    
    project_tri(dt);
    
    // Support this for T-junction glue.
    VectorXd Th_tri=Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(triangle_thickness.data(), triangle_thickness.size());
    
    advect_triangle_viaFlux( dt, Th_tri );
    clamp(Th_tri,min_th_nm,max_th_nm);
    
    triangle_thickness=std::vector<double>(Th_tri.data(),Th_tri.data()+Th_tri.size());
    
    advect_triangleVelocity( dt, triangleAdvectType);
    project_tri(dt);
    
    // Migrate triangle thickness values to vertices
    auto vertex2TriangleNaive=[this](const VectorXd& vertexValues, VectorXd& triangleValues){
        scalarVert2tri( vertexValues, triangleValues );
    };
    
    VectorXd Th_vert;
    
    auto triangle2VertexNaive=[this,num_v](const VectorXd& triangleValues, VectorXd& vertexValues){
        // take weighted average of incident triangles
        
        vertexValues.setZero(num_v);
        for(size_t vi=0;vi<num_v;++vi){
            //vertexValues[vi]=;
            const auto& inc_tris=mesh().m_vertex_to_triangle_map[vi];
            double triple_vertex_area=0.0;
            for(const auto& inc_tri:inc_tris){
                vertexValues[vi]+=triangle_areas[inc_tri]*triangleValues[inc_tri];
                triple_vertex_area+=triangle_areas[inc_tri];
            }
            vertexValues[vi]/=triple_vertex_area;
        }
    };
    triangle2VertexNaive(Th_tri,Th_vert);
    setThv(Th_vert);
    
    return dt;
}

double HGF::step_FES( double dt ) {
    
    switch(equation){
        case EulerEquation:
            actual_dt=step_EulerEquation(dt);
            break;
        case SoapFilmEquations:
            actual_dt=step_SoapFilmEquations(dt);
            break;
    }
    
    return actual_dt;
}

void HGF::init_FES_preserveVelocities(){
    Eigen::MatrixXd TriangleVelocities=ftmat();
    Eigen::MatrixXd VertexVelocities=fvmat();
    init_FES();
    setft(TriangleVelocities);
    setfv(VertexVelocities);
    setTriVelocity(TriangleVelocities);
}

void HGF::init_FES() {
    
    std::cout<<"INIT FES"<<std::endl;
    
    using namespace Eigen;
    size_t num_t = mesh().nt();
    size_t num_e = mesh().ne();
    size_t num_v = mesh().nv();
        
    triangle_centers = std::vector<Eigen::Vector3d>( num_t, Eigen::Vector3d( 0, 0, 0 ) );
    triangle_normals = std::vector<Eigen::Vector3d>( num_t, Eigen::Vector3d( 0, 0, 0 ) );
    triangle_areas = std::vector<double>( num_t, 0.0 );
    
    
    triangle_velocities=std::vector<Eigen::Vector3d>( num_t, Eigen::Vector3d( 0, 0, 0 ) );
    
    
    triangle_thickness = std::vector<double>( num_t, 0.0 );
    
    edge_centers = std::vector<Eigen::Vector3d>( num_e, Eigen::Vector3d( 0, 0, 0 ) );
    edge_t_normals = std::vector<Eigen::Vector3d>( num_e, Eigen::Vector3d( 0, 0, 0 ) );
    edge_v_normals = std::vector<Eigen::Vector3d>( num_e, Eigen::Vector3d( 0, 0, 0 ) );
    
    edge_lengths = std::vector<double>( num_e, 0.0 );
    edge_tnormal_velocities = std::vector<double>( num_e, 0.0 );
    edge_thickness = std::vector<double>( num_e, 0.0 );
    
    vertex_normals = std::vector<Eigen::Vector3d>( num_v, Eigen::Vector3d( 0, 0, 0 ) );
    
    //vertex_normals
    for (size_t vi = 0; vi < num_v; ++vi) {
        vertex_normals[ vi ] = vc( surfTrack()->get_vertex_normal( vi ) );
    }
    
    //triangle centers
    for (size_t ti = 0; ti < num_t; ++ti) {
        
        triangle_centers[ ti ] = vc( surfTrack()->get_triangle_barycenter( ti ) );
        
    }
    
    //triangle normals
    for (size_t ti = 0; ti < num_t; ++ti) {
        
        triangle_normals[ ti ] = vc( surfTrack()->get_triangle_normal( ti ) );
    }
    
    //edge_centers
    for (size_t e = 0; e < surfTrack()->m_mesh.ne(); ++e) {
        size_t v0 = mesh().m_edges[ e ][ 0 ];
        size_t v1 = mesh().m_edges[ e ][ 1 ];
        
        const LosTopos::Vec3d & x0 = surfTrack()->pm_positions[ v0 ];
        const LosTopos::Vec3d & x1 = surfTrack()->pm_positions[ v1 ];
        
        edge_centers[ e ] = vc( x0 + x1 ) / 2.0;
    }
    
    //triangle_areas
    for (size_t ti = 0; ti < num_t; ++ti) {
        
        triangle_areas[ ti ] = surfTrack()->get_triangle_area( ti );
        
    }
    
    //compute edge lengths
    for (size_t ei = 0; ei < num_e; ++ei) {
        
        edge_lengths[ ei ] = surfTrack()->get_edge_length( ei );
        
    }

    //compute edge normals
    //A. take average of the normals of two incident vertices,
    //then take cross product with edge_vertical_normal and edge.
    //B. connect centers of two incident triangles.
    const bool use_triangleNormal_or_vertexNormal = 0;
    enum EdgeNormalType{
        FromTriangleNormals,
        FromVertexNormals,
    } edgeNormalType=FromVertexNormals;
    
    for (size_t ei = 0; ei < num_e; ++ei) {
        
        size_t vi0 = mesh().m_edges[ ei ][ 0 ];
        size_t vi1 = mesh().m_edges[ ei ][ 1 ];
        
        const auto& v0 = surfTrack()->get_position( vi0 );
        const auto& v1 = surfTrack()->get_position( vi1 );
        
        const Vector3d edge = vc( v1 - v0 );
        auto&edge_normal = edge_t_normals[ ei ];
        
        switch (edgeNormalType) {
            case FromTriangleNormals:{
                size_t ti0 = mesh().m_edge_to_triangle_map[ ei ][ 0 ];
                size_t ti1 = mesh().m_edge_to_triangle_map[ ei ][ 1 ];
                
                //edge_normal = (triangle_normals[ ti1 ] + triangle_normals[ ti0 ]).normalized().cross(edge);
                
                edge_normal = (triangle_normals[ ti1 ] - triangle_normals[ ti0 ]);
                
                edge_normal.normalize();
                break;
            }
            case FromVertexNormals:{
                const Vector3d& vn0 = vertex_normals[ vi0 ];
                const Vector3d& vn1 = vertex_normals[ vi1 ];
                
                Vector3d e_vertiacal_normal = ( vn0 + vn1 ) / 2.0;
                e_vertiacal_normal -= e_vertiacal_normal.dot( edge ) / edge.squaredNorm()*edge;
                
                edge_normal = edge.cross( e_vertiacal_normal );
                
                edge_normal.normalize();
                break;
            }
        }

        //edge_v_normals
        edge_v_normals[ ei ] = edge.cross( edge_t_normals[ ei ] );
        edge_v_normals[ ei ].normalize();
        size_t adj_tri0 = mesh().m_edge_to_triangle_map[ ei ][ 0 ];
        if (edge_v_normals[ ei ].dot( triangle_normals[ adj_tri0 ] ) < 0) {
            edge_v_normals[ ei ] *= -1;
        }
        
        
    }
    
    compute_vertex_areas( vertex_areas );
    
    // Compute GradB
    GradB.setZero( 3 * num_t, 3 );
    for( size_t ti = 0; ti < num_t; ++ti )
    {
        const LosTopos::Vec3st& verts = mesh().m_tris[ ti ];
        
        for(size_t vi_t=0;vi_t<3;++vi_t)
        {
            
            const size_t& v_a = verts[ ( vi_t+1 )%3 ];
            const size_t& v_b = verts[ ( vi_t+2 ) % 3 ];
            
            const double lengthOppositeEdge = (pos(v_a)-pos(v_b )).norm();
            
            const double mag = lengthOppositeEdge/(2.0*triangle_areas[ti]);
            
            Vector3d dir=(triangle_normals[ti].cross(pos(v_b)-pos(v_a))).normalized();
            
            //flip if dir is opposite direction.
            if(dir.dot(edge_centers[mesh().get_edge_index(v_a, v_b)]-triangle_centers[ti])>0){
                dir*=-1.0;
            }
            
            GradB.row( 3 * ti + vi_t ) =mag * dir;
        }
    }
    
    DivV_vert = VectorXd::Zero( num_v );
    DivTV_vert = VectorXd::Zero( num_v );
    GradT_vert=MatrixXd::Zero( num_v,3 );
    TestVector_vert=VectorXd::Zero( num_v );
    
    std::cout<<"Compute one ring neighbors.\n"<<std::endl;
    
    mapVerticesIrregular.clear();
    one_ring_neighbors.clear();
    one_ring_neighbors.resize(num_v);
    for (size_t vi = 0; vi < num_v; ++vi) {
        if(mesh().is_vertex_nonmanifold(vi) or surfTrack()->vertex_is_any_solid(vi)){continue;}
        compute_one_ring_vertices(vi, one_ring_neighbors[vi]);
    }
    
    std::cout<<"Compute Laplacian"<<std::endl;
    // Prepare sparse solver
    computeWeightedLaplacian(MeshLaplacian);
    
    // Register ghost and real (non-ghost) triangles.
    real_triangles.clear();
    real_triangles.reserve( num_t );
    
    ghost_triangles.clear();
    ghost_triangles.reserve( num_t );
    triangles_to_real=std::vector<int>(num_t,-1);
    ghost_triangles.reserve( num_t );
    
    for( size_t ti = 0; ti < num_t; ++ti ){
        if( surfTrack()->triangle_is_all_solid( ti ) ){
            ghost_triangles.push_back( ti );
        }else{
            triangles_to_real[ti]=real_triangles.size();
            real_triangles.push_back( ti );
        }
    }
    
    // Register ghost, real, non-solid vertices
    real_vertices.clear();
    real_vertices.reserve(num_v);
    verts_to_real=std::vector<int>(num_v,-1);
    ghost_vertices.clear();
    ghost_vertices.reserve(num_v);
    
    non_solid_vertices.clear();
    non_solid_vertices.reserve(num_v);
    verts_to_non_solids=std::vector<int>(num_v,-1);
    for(size_t vi=0;vi<num_v;++vi){
        std::vector<size_t> adj_verts;
        mesh().get_adjacent_vertices(vi, adj_verts);
        
        const auto& inc_tris=mesh().m_vertex_to_triangle_map[vi];
        
        bool is_ghost_vert=true;
        for(size_t inc_tri:inc_tris){
            if(not surfTrack()->triangle_is_all_solid(inc_tri)){
                is_ghost_vert=false;
                break;
            }
        }
        
        if(is_ghost_vert){
            ghost_vertices.push_back(vi);
        }else{
            verts_to_real[vi]=real_vertices.size();
            real_vertices.push_back(vi);
        }
        
        if(not surfTrack()->vertex_is_any_solid(vi)){
            verts_to_non_solids[vi]=non_solid_vertices.size();
            non_solid_vertices.push_back(vi);
        }
        
    }
    
    // Make V_in and F_in
    computeWeightedLaplacian(MeshLaplacianIn, Eigen::VectorXd(),true);
    
    std::cout << "Compute Laplacian In" << std::endl;
    // Prepare sparse solver
    
    MeshLaplacianTransposeIn = MeshLaplacianIn.transpose();
    MLTMLIn = MeshLaplacianTransposeIn * MeshLaplacianIn;
    
    LaplacianSolverIn.analyzePattern( MLTMLIn );
    LaplacianSolverIn.factorize( MLTMLIn );
        
    // Record liquid volumes of connected components.
    numConnectedComponents=findConnectedComponents(regionsToComponents);
    std::cout<<"Num of connected components:"<<numConnectedComponents<<"\n";
    for(const auto& regionAndComponent:regionsToComponents){
        std::cout<<"region :"<<regionAndComponent.first<<" component:"<<regionAndComponent.second<<"\n";
    }
    
    initialLiquidVolumes.setZero(numConnectedComponents);
    computeLiquidVolumes(initialLiquidVolumes);
    for(const auto& regionAndComponent:regionsToComponents){
        const double volume = regionAndComponent.second!=-1?initialLiquidVolumes[regionAndComponent.second]:0;
        std::cout<<"region :"<<regionAndComponent.first<<" volume:"<<volume<<"\n";
    }
    
    std::cout<<"total liquid volume:"<<computeTotalLiquidVolume()<<"\n";
    
    setUnrenderedTrianglesAndEdges();
}

void HGF::setUnrenderedTrianglesAndEdges(){
    const size_t num_t=mesh().nt();
    const size_t num_e=mesh().ne();
    unrendered_triangles.clear();
    for(size_t ti=0;ti<num_t;++ti){
        const auto& labels=mesh().get_triangle_label(ti);
        if(labels[0]==unrendered_region and labels[1]==AIR){
            unrendered_triangles[ti]=true;
        }
    }
    
    unrendered_edges.clear();
    for(size_t ei=0;ei<num_e;++ei){
        const auto& inc_tris=mesh().m_edge_to_triangle_map[ei];
        
        if(inc_tris.size()!=2){continue;}
        
        const auto& labels0=mesh().get_triangle_label(inc_tris[0]);
        const auto& labels1=mesh().get_triangle_label(inc_tris[1]);
        
        if( labels0[0]==unrendered_region and labels0[1]==AIR
           and
           labels1[0]==unrendered_region and labels1[1]==AIR
           ){
            unrendered_edges[ei]=true;
        }
        
    }
}

void HGF::thEdge2Triangle() {
    
    int num_t = mesh().nt();
    int num_e = mesh().ne();
    //Interpolate triangle velocities from edge velocities.
    for (int tj = 0; tj < num_t; ++tj) {
        const auto& inc_edges = mesh().m_triangle_to_edge_map[ tj ];
        
        triangle_thickness[ tj ] = ( edge_thickness[ inc_edges[ 0 ] ] + edge_thickness[ inc_edges[ 1 ] ] + edge_thickness[ inc_edges[ 2 ] ] ) / 3.0;
    }
}

void HGF::velocityVertex2Triangle(bool local_energy_preserve) {
    enum VelocityVertex2TriangleType{
        Standard,
        BFECC
    }velocityVertex2TriangleType=BFECC;
    switch (velocityVertex2TriangleType) {
        case Standard:
            velocityVertex2TriangleStandard(local_energy_preserve);
            break;
        case BFECC:
            velocityVertex2TriangleBFECC(local_energy_preserve);
            break;
    }
}

void HGF::velocityTriangle2Vertex( bool local_energy_preserve ) {
    enum VelocityTriangle2VertexType{
        Standard,
        BFECC
    }velocityTriangle2VertexType=Standard;
    switch (velocityTriangle2VertexType) {
        case Standard:
            velocityTriangle2VertexStandard(local_energy_preserve);
            break;
        case BFECC:
            // For unknow reason, using BFECC with local_energy_preserve=true leads to a strange behavior.
            velocityTriangle2VertexBFECC(false);
            break;
    }
}

void HGF::velocityTriangle2VertexBFECC( bool local_energy_preserve ) {
    using namespace Eigen;
    
    // Save original triangle velocities
    MatrixXd OldVelocities=getTriVelocity();
    velocityTriangle2VertexStandard(local_energy_preserve);
    velocityVertex2TriangleStandard(local_energy_preserve);
    MatrixXd CurrentVelocities=getTriVelocity();
    
    setTriVelocity(OldVelocities+0.5*(OldVelocities-CurrentVelocities));
    velocityTriangle2VertexStandard(local_energy_preserve);
    
}

void HGF::velocityTriangle2VertexStandard( bool local_energy_preserve ) {
    
    Eigen::MatrixXd vertVelocity;
    vectorTri2vert(getTriVelocity(), vertVelocity,local_energy_preserve);
    setfv(vertVelocity);
}

// Reference: [Shi and Yu] Figure6
void HGF::interpolateVelocityFromTriangleShiYu(size_t tri_ind , const Eigen::Vector3d&position, const Eigen::MatrixXd& sourceMat, Eigen::Vector3d&velocity ) {
    
    using namespace Eigen;
    double eps = 1.0e-8;
    
    const Vector3d c0 = triangle_centers[tri_ind];
    
    // Find v0, v2 and T1.
    int v0_ind=-1, v2_ind=-1, adj_tri_ind=-1;
    const auto& inc_edges = mesh().m_triangle_to_edge_map[ tri_ind ];
    
    for (int ei_t = 0; ei_t < 3; ei_t++) {
        size_t edge = inc_edges[ ei_t ];
        
        size_t vi0 = this->mesh().m_edges[ edge ][ 0 ];
        size_t vi1 = this->mesh().m_edges[ edge ][ 1 ];
        
        const Vector3d v0 = vc(surfTrack()->get_position( vi0 ));
        const Vector3d v1 = vc(surfTrack()->get_position( vi1 ));
        
        Eigen::Vector3d hit_point;
        double u, v;
        const double EPSILON=1e-4;
        if(rayIntersectsTriangle(v0, v1, c0, position, hit_point, u, v, EPSILON) > 0){
            v0_ind = vi0;
            v2_ind = vi1;
            adj_tri_ind = mesh().m_edge_to_triangle_map[ edge ][ 0 ] != tri_ind ? mesh().m_edge_to_triangle_map[ edge ][ 0 ] : mesh().m_edge_to_triangle_map[ edge ][ 1 ];
            break;
        }
    }
    
    if(v0_ind==-1 or v2_ind==-1 or adj_tri_ind==-1){
        std::cout<<"Invalid index in interpolateVelocityFromTriangleShiYu."<<std::endl;
        std::cout<<"Consider using relaxed EPSILON for rayIntersectsTriangle e.g. 1e-3."<<std::endl;
        std::cout<<"Otherwise use \'while\' scheme for EPSILON values or just give some v0_ind, v2_ind, adj_tri_ind."<<std::endl;
        exit(-1);
    }
    
    const Vector3d v0 = vc( surfTrack()->get_position( v0_ind ) );
    const Vector3d v2 = vc( surfTrack()->get_position( v2_ind ) );
    const Vector3d c1 = triangle_centers[adj_tri_ind];
    Quaterniond rotation = RotationBetweenVectors( triangle_normals[ adj_tri_ind ], triangle_normals[ tri_ind ] );
    
    std::unordered_map <size_t,Eigen::Vector3d> modified_ad_positions;
    // v0 velicity on T0
    Vector3d vel_v0_t0 = Vector3d::Zero();
    calcVertexVelocity(v0_ind, tri_ind, vel_v0_t0,sourceMat,modified_ad_positions);
    
    // v0 velicity on T1
    // rotate velocity v0 so that it is on T0
    Vector3d vel_v0_t1 = Vector3d::Zero();
    calcVertexVelocity(v0_ind, adj_tri_ind, vel_v0_t1,sourceMat,modified_ad_positions);
    Vector3d rotated_vel_v0_t1 = rotation*vel_v0_t1;
    
    // v2 velicity on T0
    Vector3d vel_v2_t0 = Vector3d::Zero();
    calcVertexVelocity(v2_ind, tri_ind, vel_v2_t0,sourceMat,modified_ad_positions);
    
    // v2 velicity on T1
    // rotate velocity v2 so that it is on T0
    Vector3d vel_v2_t1 = Vector3d::Zero();
    calcVertexVelocity(v2_ind, adj_tri_ind, vel_v2_t1,sourceMat,modified_ad_positions);
    Vector3d rotated_vel_v2_t1 = rotation*vel_v2_t1;
    
    // Compute v0 velocity
    double small_theta_v0 = acos(clamp(((c0 - v0).normalized()).dot((position - v0).normalized()), -1.0+eps, 1.0-eps));
    double large_theta_v0 = acos(clamp(((c0 - v0).normalized()).dot((v2 - v0).normalized()), -1.0+eps, 1.0-eps))
    + acos(clamp(((v2 - v0).normalized()).dot((c1 - v0).normalized()), -1.0+eps, 1.0-eps));
    assert(std::fabs(large_theta_v0) > 1.0e-6);
    double ratio_v0 = small_theta_v0 / large_theta_v0;
    Vector3d vel_v0 = (1.0 - ratio_v0)*vel_v0_t0 + ratio_v0 * rotated_vel_v0_t1;
    
    // Compute v2 velocity
    double small_theta_v2 = acos(clamp(((c0 - v2).normalized()).dot((position - v2).normalized()), -1.0+eps, 1.0-eps));
    double large_theta_v2 = acos(clamp(((c0 - v2).normalized()).dot((v0 - v2).normalized()), -1.0+eps, 1.0-eps))
    + acos(clamp(((v0 - v2).normalized()).dot((c1 - v2).normalized()), -1.0+eps, 1.0-eps));
    assert(std::fabs(large_theta_v2) > 1.0e-6);
    double ratio_v2 = small_theta_v2 / large_theta_v2;
    Vector3d vel_v2 = (1.0 - ratio_v2)*vel_v2_t0 + ratio_v2 * rotated_vel_v2_t1;
    
    // compute p velocity
    std::vector<double> weights(3);
    
    WeightType weightType=UV;
    switch(weightType){
        case UV:{
            Vector3d dummyHitpoint;
            rayIntersectsTriangle(c0,v0,v2, position, dummyHitpoint, weights[0], weights[1]);
            weights[2]=1.0-weights[0]-weights[1];
            break;
        }
        case BaryCentric:{
            const double  t0_area = 0.5 * ( ( v0 - c0 ).cross( v2 - c0 ) ).norm();
            weights[0]= 0.5 * ( ( c0 - position ).cross( v2 - position ) ).norm()/t0_area;
            weights[1] = 0.5 * ( ( v0 - position ).cross( c0 - position ) ).norm()/t0_area;
            weights[2] = 0.5 * ( ( v2 - position ).cross( v0 - position ) ).norm()/t0_area;
            break;
        }
    }
    
    assert(weights[0]>-eps);
    assert(weights[1]>-eps);
    assert(weights[2]>-eps);
    
    velocity = (weights[0]*vel_v0 + weights[1]*vel_v2 + weights[2]*triangle_velocities[tri_ind] ) ;
    
    projectToTangent(velocity, triangle_normals[tri_ind], true);
    
}

void HGF::InterpolateScalarFromVerticesOfMultipleTriangles(const std::vector<double>& weights, const std::vector<Eigen::Vector3d>& positions, const std::vector<size_t>& tri_inds, const Eigen::VectorXd& sourceVec, double& destValue){
    
    const size_t num_points = weights.size();
    destValue=0.0;
    for(size_t p=0;p<num_points;++p){
        double partialDestValue=0.0;
        interpolateScalarFromVertices(tri_inds[p],positions[p],sourceVec,partialDestValue);
        destValue+=weights[p]*partialDestValue;
    }
}

void HGF::interpolateScalarFromVertices(size_t tri_ind , const Eigen::Vector3d&position, const Eigen::VectorXd& sourceVec, double& destValue){
    
    using namespace Eigen;
    
    int vertex0_ind = mesh().m_tris[ tri_ind ][ 0 ];
    int vertex1_ind = mesh().m_tris[ tri_ind ][ 1 ];
    int vertex2_ind = mesh().m_tris[ tri_ind ][ 2 ];
    
    std::array<double,3> weights;
    vertexWeights(weights, tri_ind, position, UV);
    
    destValue = ( weights[0]*sourceVec[vertex0_ind] + weights[1]*sourceVec[vertex1_ind]  + weights[2]*sourceVec[vertex2_ind]  ) ;
    
}

// performs the actual interpolation of the velocity
// triangle_indices, previous_positions, weights come from backtracing, local_weights are coordinates within triangles,
// sourceMat is num_t x 3 matrix that stores triangle velocities
// destVelocity is the output value, local_energy_preserve flags an energy preservation method of interpolation
void HGF::interpolateVectorFromVertices(size_t tri_ind , const std::array<double,3>& local_weights, const Eigen::MatrixXd& sourceMat, Eigen::Vector3d& destVelocity, bool local_energy_preserve){
    
    using namespace Eigen;
    
    const size_t vertex0_ind = mesh().m_tris[ tri_ind ][ 0 ];
    const size_t vertex1_ind = mesh().m_tris[ tri_ind ][ 1 ];
    const size_t vertex2_ind = mesh().m_tris[ tri_ind ][ 2 ];
    enum InterpolationType{
        naive,
        geodesicPolarMaps
    }interpolationType=geodesicPolarMaps;
    
    switch (interpolationType) {
        case naive:{
            destVelocity = ( local_weights[0]* sourceMat.row( vertex0_ind ) + local_weights[1]* sourceMat.row( vertex1_ind ) + local_weights[2]* sourceMat.row( vertex2_ind ) ) / triangle_areas[ tri_ind ];
            break;
        }
        case geodesicPolarMaps:{
            const bool simple_normal = 0;
            Vector3d pos_normal;
            if (!simple_normal) {
                pos_normal << ( local_weights[0]*vertex_normals[ vertex0_ind ] + local_weights[1]*vertex_normals[ vertex1_ind ] + local_weights[2]*vertex_normals[ vertex2_ind ] );
                pos_normal.normalize();
            }
            else {
                pos_normal << triangle_normals[ tri_ind ];
            }
            //project ver_vel0,1,2 to tangent plane
            Vector3d projected_vel0=sourceMat.row( vertex0_ind );
            if(mesh().is_vertex_nonmanifold(vertex0_ind)){
                projectToTangent(projected_vel0, vertex_normals[vertex0_ind], pos_normal);
            }
            Vector3d projected_vel1=sourceMat.row( vertex1_ind );
            if(mesh().is_vertex_nonmanifold(vertex1_ind)){
                projectToTangent(projected_vel1, vertex_normals[vertex1_ind], pos_normal);
            }
            Vector3d projected_vel2=sourceMat.row( vertex2_ind );
            if(mesh().is_vertex_nonmanifold(vertex2_ind)){
                projectToTangent(projected_vel1, vertex_normals[vertex2_ind], pos_normal);
            }
            
            destVelocity = local_weights[0]*projected_vel0 + local_weights[1]*projected_vel1 + local_weights[2]*projected_vel2 ;
            break;
        }
    }
    
    if ( local_energy_preserve) {
        double local_energy = 0.0;
        local_energy += local_weights[0]* sourceMat.row( vertex0_ind ).squaredNorm();
        local_energy += local_weights[1]* sourceMat.row( vertex1_ind ).squaredNorm();
        local_energy += local_weights[2]* sourceMat.row( vertex2_ind ).squaredNorm();
        
        if(!destVelocity.isZero()){
            destVelocity.normalize();
            destVelocity *= std::pow( local_energy, 0.5 );
        }
    }
}

// triangle_indices, previous_positions, weights come from backtracing, sourceMat is num_t x 3 matrix that stores triangle velocities
// destVelocity is the output value, local_energy_preserve flags an energy preservation method of interpolation
void HGF::interpolateVectorFromVertices(size_t tri_ind , const Eigen::Vector3d& position, const Eigen::MatrixXd& sourceMat, Eigen::Vector3d& destVelocity, bool local_energy_preserve){
    
    std::array<double,3> local_weights;
    vertexWeights(local_weights, tri_ind, position, UV);
    
    interpolateVectorFromVertices(tri_ind, local_weights, sourceMat, destVelocity, local_energy_preserve);
}

// triangle_indices, previous_positions, weights come from backtracing, sourceMat is num_t x 3 matrix that stores vertex velocities
// destVelocity is the output value, local_energy_preserve flags an energy preservation method of interpolation
void HGF::interpolateVectorFromVerticesOnWeightedPaths(const std::vector<size_t>& tri_ind , const std::vector<Eigen::Vector3d>& position, const std::vector<double>& backtrace_weights, const Eigen::Matrix3d& begin_local_frame, const std::vector<Eigen::Matrix3d>& end_local_frames, const Eigen::MatrixXd& sourceMat, Eigen::Vector3d& destVelocity, bool local_energy_preserve) {
    
    const size_t num_triangles = backtrace_weights.size();
    std::array<double,3> local_weights;
    
    assert((begin_local_frame != Eigen::Matrix3d::Zero()) && (begin_local_frame != Eigen::Matrix3d::Zero()));
    
    for(size_t p=0;p<num_triangles;++p){
        Eigen::Vector3d partialDestVelocity( 0.0, 0.0, 0.0 );
        vertexWeights(local_weights, tri_ind[p], position[p], UV);
        interpolateVectorFromVertices(tri_ind[p], local_weights, sourceMat, partialDestVelocity, local_energy_preserve);
        
        const double cx = partialDestVelocity.dot(end_local_frames[p].row(0));
        const double cy = partialDestVelocity.dot(end_local_frames[p].row(1));
        const double cz = partialDestVelocity.dot(end_local_frames[p].row(2));
        
        partialDestVelocity = cx * begin_local_frame.row(0) + cy * begin_local_frame.row(1) + cz * begin_local_frame.row(2);
        
        destVelocity += backtrace_weights[p]*partialDestVelocity;
    }
}

// Given a triangle index and a position, figure out the local coordinate that define the given position within the given triangle
// return the vector of weights
void HGF::vertexWeights(std::array<double,3>& weights, size_t tri_ind, const Eigen::Vector3d&position, WeightType weightType){
    using namespace Eigen;
    
    switch(weightType){
        case UV:{
            const auto& tri=mesh().m_tris[tri_ind];
            Vector3d dummyHitpoint;
            double u,v,w;
            // dummyHitPoint, u,v are output values
            rayIntersectsTriangle(pos(tri[0]), pos(tri[1]), pos(tri[2]), position, dummyHitpoint, u, v);
            w=1.0-u-v;
            weights[0]=w;
            weights[1]=u;
            weights[2]=v;
            break;
        }
        case BaryCentric:{
            int vertex0_ind = mesh().m_tris[ tri_ind ][ 0 ];
            int vertex1_ind = mesh().m_tris[ tri_ind ][ 1 ];
            int vertex2_ind = mesh().m_tris[ tri_ind ][ 2 ];
            
            const Vector3d vertex0 = vc( surfTrack()->get_position( vertex0_ind ) );
            const Vector3d vertex1 = vc( surfTrack()->get_position( vertex1_ind ) );
            const Vector3d vertex2 = vc( surfTrack()->get_position( vertex2_ind ) );
            
            weights[0] = 0.5 * ( ( vertex2 - position ).cross( vertex1 - position ) ).norm()/triangle_areas[tri_ind];
            weights[1] = 0.5 * ( ( vertex0 - position ).cross( vertex2 - position ) ).norm()/triangle_areas[tri_ind];
            weights[2] = 0.5 * ( ( vertex1 - position ).cross( vertex0 - position ) ).norm()/triangle_areas[tri_ind];
            break;
        }
    }
}

void HGF::thTriangle2Edge() {
    using namespace Eigen;
    int num_e = mesh().ne();
    
    const bool simple_averaging = 0;
    
    for (int ei = 0; ei < num_e; ++ei) {
        //
        auto inc_tris = mesh().m_edge_to_triangle_map[ ei ];
        size_t t0 = inc_tris[ 0 ];
        size_t t1 = inc_tris[ 1 ];
        
        if (simple_averaging) {
            //Averaging velocities of incident triangles.
            edge_thickness[ ei ] = triangle_thickness[ t0 ] + triangle_thickness[ t1 ] / 2.0;
        }
        else {
            //Averaging based on the distances between edge and triangles.
            double dist0 = ( triangle_centers[ t0 ] - edge_centers[ ei ] ).norm();
            double dist1 = ( triangle_centers[ t1 ] - edge_centers[ ei ] ).norm();
            edge_thickness[ ei ] = ( dist1*triangle_thickness[ t0 ] + dist0*triangle_thickness[ t1 ] ) / ( dist0 + dist1 );
        }
        
    }
}

void HGF::velocityVertex2TriangleStandard(bool local_energy_preserve) {
    Eigen::MatrixXd triMat;
    vectorVert2tri(fvmat(), triMat, local_energy_preserve);
    setTriVelocity(triMat);
}

void HGF::velocityVertex2TriangleBFECC(bool local_energy_preserve) {
    using namespace Eigen;
    // Save original vertex velocities
    MatrixXd OldVelocities=fvmat();
    velocityVertex2TriangleStandard(local_energy_preserve);
    velocityTriangle2VertexStandard(local_energy_preserve);
    MatrixXd CurrentVelocities=fvmat();
    setfv(OldVelocities+0.5*(OldVelocities-CurrentVelocities));
    velocityVertex2TriangleStandard(local_energy_preserve);
    
}

void HGF::velocityTriangle2VertexFLIP(Eigen::MatrixXd& OldTriVelocities, Eigen::MatrixXd&OldVertVelocities){
    
    using namespace Eigen;
    const size_t num_t=mesh().nt();
    const size_t num_v=mesh().nv();
    
    // Set delta_triVel to get delta_vertVel
    Eigen::MatrixXd DeltaTriVelocities=getTriVelocity()-OldTriVelocities;
    for(size_t ti=0;ti<num_t;++ti){
        Vector3d DeltaVelocity=DeltaTriVelocities.row(ti);
        projectToTangent(DeltaVelocity, triangle_normals[ti]);
        DeltaTriVelocities.row(ti)=DeltaVelocity;
    }
    setTriVelocity(DeltaTriVelocities);
    
    // Convert DeltaTriVelocities to OldVertVelocities
    velocityTriangle2Vertex(true);
    // Add DeltaVertVelocities to fv.
    MatrixXd DeltaVertVelocities = fvmat();
    for(size_t vi=0;vi<num_v;++vi){
        Vector3d OldVelocity=OldVertVelocities.row(vi);
        projectToTangent(OldVelocity, vertex_normals[vi]);
        OldVertVelocities.row(vi)=OldVelocity;
        
        Vector3d DeltaVelocity=DeltaVertVelocities.row(vi);
        projectToTangent(DeltaVelocity, vertex_normals[vi]);
        DeltaVertVelocities.row(vi)=DeltaVelocity;
    }
    
    setfv(OldVertVelocities + DeltaVertVelocities);
}

void HGF::thVertex2Edge() {
    
    using namespace Eigen;
    
    const size_t num_e = mesh().ne();
    const size_t  num_t = mesh().nt();
    
    bool vertex_base = 1;
    if (vertex_base) {
        for (int ei = 0; ei < num_e; ++ei) {
            const auto& vertices = mesh().m_edges[ ei ];
            size_t vi0 = vertices[ 0 ];
            size_t vi1 = vertices[ 1 ];
            
            edge_thickness[ ei ] = ( thv ( vi0 ) + thv( vi1 ) ) / 2.0;
        }
    }
    else {
        Eigen::VectorXd Th_tri;
        scalarVert2tri(thvv(), Th_tri);
        
        for (int ti = 0; ti < num_t; ++ti) {
            triangle_thickness[ti]=Th_tri[ti];
        }
        
        thTriangle2Edge();
        
    }
}

void HGF::thEdge2Vertex()
{
    using namespace Eigen;
    size_t num_v = mesh().nv();
    
    bool area_aware = 1;//generally false is better
    
    if (area_aware) {
        for (size_t vi = 0; vi < num_v; ++vi) {
            const auto& edges = mesh().m_vertex_to_edge_map[ vi ];
            size_t num_ev = edges.size();
            
            double sum_th = 0.0;
            VectorXd edge_areas = VectorXd::Zero( num_ev );
            
            for (size_t ei_v = 0; ei_v < num_ev; ++ei_v) {
                size_t edge = edges[ ei_v ];
                
                // Area-aware version
                const auto & inc_tris_to_edge = mesh().m_edge_to_triangle_map[ edge ];
                size_t t0 = inc_tris_to_edge[ 0 ];
                size_t t1 = inc_tris_to_edge[ 1 ];
                edge_areas[ ei_v ] = 1 / 3.0*( triangle_areas[ t0 ] + triangle_areas[ t1 ] );
                
                sum_th += edge_areas[ ei_v ] * edge_thickness[ edge ];
                
                
            }
            
            thv( vi ) = sum_th / edge_areas.sum();
            
        }
    }
    else {
        for (size_t vi = 0; vi < num_v; ++vi) {
            const auto& edges = mesh().m_vertex_to_edge_map[ vi ];
            size_t num_ev = edges.size();
            double sum_th = 0.0;
            
            for (size_t ei_v = 0; ei_v < num_ev; ++ei_v) {
                size_t edge = edges[ ei_v ];
                
                sum_th += edge_thickness[ edge ];
                
            }
            
            thv( vi ) = sum_th / num_ev;
            
        }
    }
    
}

void HGF::computeRotationOnVertices(const Eigen::MatrixXd& sourceMat, Eigen::MatrixXd& rotMat){
    using namespace Eigen;
    
    size_t num_v = mesh().nv();
    size_t num_t = mesh().nt();
    
    rotMat.setZero(num_v,3);
    
    for( size_t ti:real_triangles )
    {
        const LosTopos::Vec3st& verts = mesh().m_tris[ ti ];
        for( size_t vi_t = 0; vi_t < 3; ++vi_t )
        {
            const size_t& v = verts[ vi_t ];

            const Vector3d GradBVec=GradB.row( ti * 3 + vi_t );
            const Vector3d SourceVec=sourceMat.row(ti) ;
            
            rotMat.row( v )+= triangle_areas[ti]*GradBVec.cross( SourceVec );
            
        }
    }
    
    for( size_t vi : real_vertices )
    {
        Vector3d rotVec=rotMat.row(vi);
        projectToTangent(rotVec, vertex_normals[vi],false);
        rotMat.row(vi)=rotVec;
    }
    
    // apply mass matrix
    MatrixXd V;
    MatrixXi F;
    meshAsMatrices(V,F);
    SparseMatrix<double> M, Minv;
    igl::massmatrix( V, F, igl::MASSMATRIX_TYPE_DEFAULT, M );
    igl::invert_diag( M, Minv );
    
    // Multiplying 0.5 to get the same result as Laplacian made via
    // libigl's cotmatrix and mass matrixv.
    rotMat=Minv*rotMat;
    
}

void HGF::computeWeightedDivergenceOnVertices(const Eigen::MatrixXd& vectorField, Eigen::VectorXd& divVec, const Eigen::VectorXd& scalarField){
    
    using namespace Eigen;
    
    size_t num_v = mesh().nv();
    size_t num_t = mesh().nt();
    
    MatrixXd sourceVector=vectorField;
    
    if(scalarField.size()==num_t){
        for(size_t ti=0;ti<num_t;++ti){
            sourceVector.row(ti)*=scalarField[ti];
        }
    }
    
    divVec.setZero( num_v );
    
    switch(divergenceMethod){
        case deGoesEtAl:{
            
            // First approach
            //            for( size_t ti = 0; ti < num_t; ++ti )
            for( size_t tti = 0; tti < real_triangles.size(); ++tti )
                
            {
                size_t ti=real_triangles[tti];
                
                const LosTopos::Vec3st& verts = mesh().m_tris[ ti ];
                for( size_t vi_t = 0; vi_t < 3; ++vi_t )
                {
                    const size_t& v = verts[ vi_t ];
                    const size_t& v_a = verts[ ( vi_t+1 )%3 ];
                    const size_t& v_b = verts[ ( vi_t+2 ) % 3 ];
                    
                    size_t e_ind=mesh().get_edge_index(v_a, v_b);
                    const Vector3d edge = pos(v_b)-pos(v_a );
                    Vector3d outerVector=triangle_normals[ti].cross(edge);
                    if(outerVector.dot(edge_centers[e_ind]-triangle_centers[ti])<0){
                        outerVector*=-1;
                    }
                    
                    divVec[ v ] +=  sourceVector.row(ti).dot(outerVector);
                }
            }
            
            divVec*=0.5;
            
            break;
        }
        case TongEtAl:{
            
            // A-2. compute Poisson.
            //
            
            // http://www.cs.kent.edu/~zwang/schedule/zy5.pdf
            // Eqn. 5
            // To facilitate the
            // implementation, we note that ∇φik is simply the vector orthogonal
            // to the face fik opposite to i in the tet Tk, pointing towards i and with
            //a magnitude of |∇φik| = area(fik) / 3 |Tk| .
            //
            // this is for tetrahedrons.
            // For triangles, we should modify
            // |∇φik| = length(fik) / (2 * |Tk|) where
            // fik is the edge opposite of vertex i and |Tk| is the area of triangle.
            // https://www.microsoft.com/en-us/research/wp-content/uploads/2016/02/poisson.pdf
            // Eqn. 4
            
            // do only real triangles
            for( size_t tti = 0; tti < real_triangles.size(); ++tti )
                
            {
                size_t ti=real_triangles[tti];
                
                const LosTopos::Vec3st& verts = mesh().m_tris[ ti ];
                for( size_t vi_t = 0; vi_t < 3; ++vi_t )
                {
                    const size_t& v = verts[ vi_t ];
                    divVec[ v ] += triangle_areas[ti]* GradB.row( ti * 3 + vi_t ).dot( sourceVector.row(ti)  ) ;
                    
                }
            }
            
            break;
        }
    }
    
    MatrixXd V;
    MatrixXi F;
    meshAsMatrices( V, F );
    SparseMatrix<double> M, Minv;
    igl::massmatrix( V, F, igl::MASSMATRIX_TYPE_DEFAULT, M );
    igl::invert_diag( M, Minv );
    
    // Multiplying 0.5 to get the same result as Laplacian made via
    // libigl's cotmatrix and mass matrixv.
    divVec = Minv * divVec;
    
}

void HGF::computePressureOnVertices(double dt, const Eigen::VectorXd& Uth, const Eigen::VectorXd& Th){
    using namespace Eigen;
    
    // Compute div V at vertices.
    const size_t num_t = mesh().nt();
    const size_t num_v = mesh().nv();
    const size_t r_numv=real_vertices.size();
    
    std::cout<<"Start Solving Pressure Matrix"<<std::endl;
    double maxError;
    
    if(Th.size()==0)
    {
        MatrixXd TriangleVelocity(num_t,3);
        for(size_t ti=0;ti<num_t;++ti){
            TriangleVelocity.row(ti)=triangle_velocities[ti];
        }
        computeWeightedDivergenceOnVertices(TriangleVelocity,DivV_vert);
        
        // Change DivV_vert to DivV_vert_in.
        DivV_vertIn.setZero(r_numv);
        for(size_t rvi = 0;rvi<r_numv;++rvi){
            DivV_vertIn[rvi]=DivV_vert[real_vertices[rvi]];
        }
        
        VectorXd Uth_in;
        if(Uth.size()>0){
            Uth_in.resize(r_numv);
            for(size_t rvi=0;rvi<r_numv;++rvi){
                Uth_in[rvi]=Uth[real_vertices[rvi]];
            }
        }
        
        const VectorXd& targetVector=Uth_in.size() > 0? DivV_vertIn+Uth_in:DivV_vertIn;
        P_vertIn=LaplacianSolverIn.solve(MeshLaplacianTransposeIn*targetVector).eval();
        maxError=(MeshLaplacianIn*P_vertIn-targetVector ).array().abs().maxCoeff();
        
        P_vert.setZero(num_v);
        for(size_t rvi = 0;rvi<real_vertices.size();++rvi){
            P_vert[real_vertices[rvi]]=P_vertIn[rvi];
        }
        
    }else{
        MatrixXd TriangleVelocity(num_t,3);
        for(size_t ti=0;ti<num_t;++ti){
            TriangleVelocity.row(ti)=triangle_velocities[ti];
        }
        
        // Thickness at triangles
        VectorXd TriangleT;
        scalarVert2tri( Th, TriangleT );
        
        VectorXd DivTv_vert;
        computeWeightedDivergenceOnVertices( TriangleVelocity, DivTv_vert, TriangleT );
        
        VectorXd DivTv_vert_in;
        DivTv_vert_in.setZero(r_numv);
        for(size_t rvi = 0;rvi<r_numv;++rvi){
            DivTv_vert_in[rvi]=DivTv_vert[real_vertices[rvi]];
        }
        VectorXd Uth_in;
        if(Uth.size()>0){
            Uth_in.resize(r_numv);
            for(size_t rvi=0;rvi<r_numv;++rvi){
                Uth_in[rvi]=Uth[real_vertices[rvi]];
            }
        }
        
        // Right-hand side of Poisson's equation
        const VectorXd& targetVector=Uth_in.size() > 0? DivTv_vert_in+Uth_in:DivTv_vert_in;
        P_vertIn=LaplacianSolverIn.solve(MeshLaplacianTransposeIn*targetVector).eval();
        maxError=(MeshLaplacianIn*P_vertIn-targetVector ).array().abs().maxCoeff();
        
        
        P_vert.setZero(num_v);
        for(size_t rvi = 0;rvi<real_vertices.size();++rvi){
            P_vert[real_vertices[rvi]]=P_vertIn[rvi];
        }
        
    }
    std::cout<<"Linear System Solving Done. Error="<<maxError<<std::endl;
    std::cout<<"Finish Solving Pressure Matrix"<<std::endl;
    
}

void HGF::th_laplacian_filter(double filter_coef){
    using namespace Eigen;
    size_t num_v=mesh().nv();
    LT_vert.setZero(num_v);
    
    enum Laplacian_Mode{
        naive,
        matrix
    };
    Laplacian_Mode laplacian_mode=matrix;
    
    switch (laplacian_mode) {
        case Laplacian_Mode::naive :{
            
            std::vector<size_t>adjacent_vertices;
            for(size_t vi=0;vi<num_v;++vi){
                mesh().get_adjacent_vertices( vi, adjacent_vertices );
                
                double area_sum = 0.0;
                double th_sum = 0.0;
                for (const auto& ad_v : adjacent_vertices) {
                    th_sum += thv( ad_v );
                }

                double diffFromNeighbors=thv( vi )-th_sum / adjacent_vertices.size();
                thv(vi)-=filter_coef*diffFromNeighbors;
                
                LT_vert[vi]=diffFromNeighbors;
                
            }
            
            break;
        }
        case Laplacian_Mode::matrix :{
            
            VectorXd Th;
            Th.setZero( num_v );
            for( size_t i = 0; i < num_v; i++ )
            {
                Th[i] = thv( i );
            }
            
            VectorXd ThFiltered =  MeshLaplacian* Th;
            LT_vert = ThFiltered;
            
            for( size_t i = 0; i < num_v; i++ )
            {
                thv( i )+=filter_coef*LT_vert[i];
            }
            
            break;
        }

    }
    
    
}

// Find the intersection of 3D line segments AB and CD with some robustness.
int HGF::twoSegmentsIntersect( Eigen::Vector3d& result, const Eigen::Vector3d& A, const Eigen::Vector3d& B, const Eigen::Vector3d& C, const Eigen::Vector3d& D ) {
    
    using namespace Eigen;
    const double epsilon = 1e-6; 1e-6; //1e-6 is standard
    
    // Fails when exactly AB = CD.
    if (( B - A ).norm() < epsilon || ( D - C ).norm() < epsilon) return 0;
    
    Vector3d AB = B - A;
    Vector3d CD = D - C;
    
    Vector3d n1 = AB.normalized();
    Vector3d n2 = CD.normalized();
    
    double work1 = n1.dot( n2 );
    double work2 = 1.0 - work1*work1;
    
    // Fails AB and CD are exactly parallel
    if (work2 < epsilon)
    {
        return 0;
    }
    
    Vector3d AC = C - A;
    
    double d1 = ( AC.dot( n1 ) - work1*AC.dot( n2 ) ) / work2;
    double d2 = ( work1*AC.dot( n1 ) - AC.dot( n2 ) ) / work2;
    
    // Nearest point on AB
    Vector3d result1 = A + d1*n1;
    //Nearest point on BC
    Vector3d result2 = C + d2*n2;
    
    // Evaluate intersection
    if (( result2 - result1 ).norm() < epsilon) {
        // When AB and CD intersect
        if (( d1 > -epsilon ) && ( d2 > -epsilon ) && ( d1 < ( ( B - A ).norm() + epsilon ) ) && ( d2 < ( ( D - C ).norm() + epsilon ) )) {
            // When AB and BC intersect
            result = result1;
            return 1;
        }
    }
    
    // When no intersection happens
    return -1;
}

// Evaluate if a given point is inside a triangle.
bool HGF::inside_triangle2d( size_t tri_ind, const Eigen::Vector3d& P, double epsilon )
{
    
    using namespace Eigen;
    
    const size_t& vi0 = mesh().m_tris[ tri_ind ][ 0 ];
    const size_t& vi1 = mesh().m_tris[ tri_ind ][ 1 ];
    const size_t& vi2 = mesh().m_tris[ tri_ind ][ 2 ];
    
    const Vector3d A = vc( surfTrack()->get_position( vi0 ) );
    const Vector3d B = vc( surfTrack()->get_position( vi1 ) );
    const Vector3d C = vc( surfTrack()->get_position( vi2 ) );
    
    Vector3d AB = B - A;
    Vector3d BP = P - B;
    
    Vector3d BC = C - B;
    Vector3d CP = P - C;
    
    Vector3d CA = A - C;
    Vector3d AP = P - A;
    
    Vector3d c1 = AB.cross( BP );
    Vector3d c2 = BC.cross( CP );
    Vector3d c3 = CA.cross( AP );
    
    double dot_12 = c1.dot( c2 );
    double dot_23 = c2.dot( c3 );
    double dot_31 = c3.dot( c1 );
    
    // Consider the case the point is on the extension of one of the triangle edges.
    if (dot_12 > -epsilon && dot_23 > -epsilon && dot_31 > -epsilon) {
        return true;
    }
    
    // The point is outside the triangle
    return false;
}

bool HGF::inside_triangle3d( size_t tri_ind, const Eigen::Vector3d& P, double EPSILON )
{
    
    using namespace Eigen;
    
    // make ray
    Vector3d hit_point;
    return rayIntersectsTriangle( tri_ind, P, triangle_normals[ tri_ind ], hit_point, EPSILON ) > 0;
    
}

int HGF::rayIntersectsTriangle( const Eigen::Vector3d& position0, const Eigen::Vector3d& position1, const Eigen::Vector3d& position2, const Eigen::Vector3d& rayOrigin, Eigen::Vector3d& hit_point, double& u, double& v, double EPSILON ) {
    
    using namespace Eigen;
    
    const Vector3d edge1 = position1 - position0;// corresponding to u
    const Vector3d edge2 = position2 - position0;// corresponding to v
    
    const Vector3d rayDirection = edge1.cross( edge2 ).normalized();    // triangle normal
    
    // h, s, q;
    const Vector3d  h = rayDirection.cross( edge2 );    // vector ortho to rayDirection & edge2 in edge1,edge2-plane
    const double a = edge1.dot( h );                    // collinearity of edge1 & h
    
    // For the case ray is parallel to the triangle, i.e. ray in on the triangle.
    if (a > -EPSILON && a < EPSILON)
    {
        hit_point = rayOrigin;
        return 1;
    }
    
    // We compute hit point (intersection) first so that the user can use no matter when this function ends.
    const double f = 1 / a;
    const Vector3d s = rayOrigin - position0;
    u = f * ( s.dot( h ) );
    
    const Vector3d q = s.cross( edge1 );
    v = f * rayDirection.dot( q );
    
    double t = f * edge2.dot( q );
    
    hit_point = rayOrigin + rayDirection * t;
    
    if (u < -EPSILON || u > 1 + EPSILON) {
        return -1;
    }
    
    if (v < -EPSILON || u + v > 1.0 + EPSILON) {
        return -2;
    }
    
    return 1;
    
}

int HGF::rayIntersectsTriangle( size_t tri_ind, const Eigen::Vector3d& rayOrigin, const Eigen::Vector3d& rayDirection, Eigen::Vector3d& hit_point, double EPSILON ) {
    
    using namespace Eigen;
    
    const size_t& vi0 = mesh().m_tris[ tri_ind ][ 0 ];
    const size_t& vi1 = mesh().m_tris[ tri_ind ][ 1 ];
    const size_t& vi2 = mesh().m_tris[ tri_ind ][ 2 ];
    
    const Vector3d vertex0 = vc( surfTrack()->get_position( vi0 ) );
    const Vector3d vertex1 = vc( surfTrack()->get_position( vi1 ) );
    const Vector3d vertex2 = vc( surfTrack()->get_position( vi2 ) );
    const Vector3d edge1 = vertex1 - vertex0;
    const Vector3d edge2 = vertex2 - vertex0;
    
    //, h, s, q;
    const Vector3d  h = rayDirection.cross( edge2 );
    const double a = edge1.dot( h );
    
    // For the case case ray is parallel to the triangle, i.e. ray in on the triangle.
    if (a > -EPSILON && a < EPSILON)
        //if(a>0 && a<0)
        //if (a > -std::numeric_limits<double>::epsilon() && a < std::numeric_limits<double>::epsilon())
    {
        hit_point = rayOrigin;
        return 1;
    }
    
    // We compute hit point (intersection) first so that the user can use no matter when this function ends.
    const double f = 1 / a;
    const Vector3d s = rayOrigin - vertex0;
    const double u = f * ( s.dot( h ) );
    
    const Vector3d q = s.cross( edge1 );
    const double v = f * rayDirection.dot( q );
    
    double t = f * edge2.dot( q );
    
    hit_point = rayOrigin + rayDirection * t;
    
    if (u < -EPSILON || u > 1 + EPSILON) {
        //    if (u < 0.0 || u > 1.0){
        return -1;
    }
    
    if (v < -EPSILON || u + v > 1.0 + EPSILON) {
        //     if (v < 0.0 || u + v > 1.0){
        return -2;
    }
    
    return 1;
    
}

// backtrace from a vertex on the first triangle (the velocity was already determined)
// input: time step, source vertex, source triangle
int HGF::backtraceFromVertex(const double dt, const size_t v_ind, const size_t tri_ind, const Eigen::Vector3d& projected_prev_velocity, const std::unordered_map <size_t,Eigen::Vector3d>& modified_ad_positions,Eigen::Vector3d& naive_prev_position,std::vector<Eigen::Vector3d>& previous_positions, std::vector<double>& weights, std::vector<size_t>& triangle_indices) {
    
// For debug purpose
#ifdef advect_experiment
    if (!intersectionOnce) {
        trace.push_back( pos( v_ind ) );
    }
#endif //advect_experiment
    

    using namespace Eigen;
    Vector3d tri_normal = triangle_normals[ tri_ind ];
    const auto& tri = mesh().m_tris[ tri_ind ];
    
    // save the two vertices from source triangle that are not the source vertex
    size_t adj_ver_ind[ 2 ];
    for (int vi = 0, i = 0; vi < 3; ++vi) {
        if (tri[ vi ] != v_ind) {
            adj_ver_ind[ i++ ] = tri[ vi ];
        }
    }
    
    const Vector3d& vert_pos = Vector3d(0, 0, 0);
    const Vector3d& projected_pos1 = modified_ad_positions.at( adj_ver_ind[ 0 ] );
    const Vector3d& projected_pos2 = modified_ad_positions.at( adj_ver_ind[ 1 ] );
    
    Vector3d hit_point;
    double u, v;
    
    // Check naive_prev_position is inside the (projected) triangle (in 2D).
    // hit_point, u,v are output values
    int inside_triangle = rayIntersectsTriangle( vert_pos, projected_pos1, projected_pos2, naive_prev_position, hit_point, u, v, 1e-6 );
    
    // For the case u~v~0, i.e. backtraced position is very close to pos(v_ind)
    // should be already caught in the beginning of semiLagrangianFromVertex function and should not reach here.
    const double eps = 1e-10;
    assert( std::abs( u ) > eps or std::abs( v ) > eps );
    
    // There are three cases.
    // 1. u>=0 && v>=0 && u+v<=1. Inside triangle.
    // 2. u>=0 && v>=0 && u+v>1. Intersects with an edge or two.
    // 3. u<0 || v<0. Does not intersect any edge.
    
    // Case 1
    if (inside_triangle == 1) {//1e-6 standard value
        
        // compute reall prev_position  using u and v
        const Vector3d pos_vert = pos( v_ind );
        const Vector3d pos1 = pos( adj_ver_ind[ 0 ] );
        const Vector3d pos2 = pos( adj_ver_ind[ 1 ] );
        previous_positions.push_back( pos_vert + u * ( pos1 - pos_vert ) + v*( pos2 - pos_vert ));
        weights.push_back(1);
        triangle_indices.push_back(tri_ind);
        
#ifdef advect_experiment
        if (!intersectionOnce) {
            trace.push_back( previous_positions[0] );
            path_lengths.push_back( ( previous_positions[0] - pos( v_ind ) ).norm() );
            std::cout << "path lengths\n";
            for (auto path : path_lengths) {
                std::cout << path << ",";
            }std::cout << "\n\n";
            std::cout << "trace\n";
            for (auto point : trace) {
                std::cout << point.x() << "," << point.y() << "," << point.z() << "\n";
            }std::cout << "\n";
            
            intersectionOnce = true;
        }
#endif //advect_experiment
        
        return 0;
    }
    // Case 3
    else if (u < 0 or v < 0) {
        return -1;
    }
    // Case 2
    else {
        
        // If we normalize u and v such that u'+v'=1,
        // then u' * (pos1-pos_vert)+ v' *(pos2-pos_vert) is a point on the edge connecting pos1 and pos2.
        const double u_prime = u / ( u + v );
        const double v_prime = v / ( u + v );
        // compute real prev_position  using u and v
        const Vector3d pos_vert = pos( v_ind );
        const Vector3d pos1 = pos( adj_ver_ind[ 0 ] );
        const Vector3d pos2 = pos( adj_ver_ind[ 1 ] );
        Vector3d segsIntersection = pos_vert +u_prime * (pos1- pos_vert) + v_prime*(pos2- pos_vert);
        
        const auto& inc_edges = mesh().m_triangle_to_edge_map[ tri_ind ];
        // edge opposite of the vertex v_ind
        
        int opposite_edge_ind = -1;
        for (size_t ei_t = 0; ei_t < 3; ++ei_t) {
            const auto& edge_ind = inc_edges[ ei_t ];
            if (mesh().m_edges[ edge_ind ][ 0 ] != v_ind and mesh().m_edges[ edge_ind ][ 1 ] != v_ind) {
                opposite_edge_ind = edge_ind;
            }
        }
        
        assert( opposite_edge_ind != -1 );
        
        if( surfTrack()->edge_is_all_solid( opposite_edge_ind ) )
        {
#ifdef advect_experiment
            if( !intersectionOnce ) {
                trace.push_back( segsIntersection );
                intersectionOnce = true;
            }
#endif //advect_experiment
            previous_positions.push_back( segsIntersection );
            weights.push_back(1);
            triangle_indices.push_back(tri_ind);
            return 0;
        }
        
        int backtrace_success=-1;
        
        const Vector3d projected_prev_position=pos_vert-dt*projected_prev_velocity;
        
        // If we are on a T-junction, decide if we should split path for each adjacent triangles.
        const bool is_on_Tjunction=mesh().is_edge_nonmanifold(opposite_edge_ind);
        std::vector<size_t> valid_adj_triangles;
        std::vector<double> incoming_water_rates;
        if(is_on_Tjunction){
            collectValidTrianglesAndIncomingWaters(dt,tri_ind,opposite_edge_ind,valid_adj_triangles, incoming_water_rates);
            
            // Finish backtrace if there are no triangles from which water is incoming.
            if(valid_adj_triangles.empty()){
                previous_positions.push_back(segsIntersection);
                weights.push_back(1);
                triangle_indices.push_back(tri_ind);
                
#ifdef advect_experiment
                if( !intersectionOnce ) {
                    trace.push_back( segsIntersection );
                    intersectionOnce = true;
                }
#endif //advect_experiment
                return -2;
            }
            
            for(size_t ti_adj=0;ti_adj<valid_adj_triangles.size();++ti_adj){
                const size_t adj_tri_ind=valid_adj_triangles[ti_adj];
                
                const Vector3d backtraceDirection=(segsIntersection-pos_vert).normalized();
                const double backtraceLength=dt*projected_prev_velocity.norm();
                Vector3d next_projected_prev_position=
                findNextProjectedPrevPosition(tri_ind,adj_tri_ind,pos_vert +backtraceLength*backtraceDirection,segsIntersection,opposite_edge_ind);
                
                std::vector<Vector3d> partial_previous_positions;
                std::vector<double> partial_weights;
                std::vector<size_t> partial_triangle_indices;
                
                backtrace_success = backtraceFromEdge_woFrame(dt, opposite_edge_ind, adj_tri_ind, segsIntersection, next_projected_prev_position, partial_previous_positions, partial_weights, partial_triangle_indices, 1);
                
                if(backtrace_success==-1){
                    
                    return -1;
                }
                
                previous_positions.insert(previous_positions.end(), partial_previous_positions.begin(),partial_previous_positions.end());
                triangle_indices.insert(triangle_indices.end(), partial_triangle_indices.begin(),partial_triangle_indices.end());
                for(const double& partial_weight: partial_weights){
                    weights.push_back(partial_weight*incoming_water_rates[ti_adj]);
                }
                
            }
            
            return -2;
        }else{
            const auto& adj_tris=mesh().m_edge_to_triangle_map[ opposite_edge_ind ];
            const size_t adj_tri_ind=adj_tris[0]!=tri_ind?adj_tris[0]:adj_tris[1];
            
            const Vector3d backtraceDirection=(segsIntersection-pos_vert).normalized();
            const double backtraceLength=std::fabs(dt)*projected_prev_velocity.norm();
            Vector3d next_projected_prev_position=
            findNextProjectedPrevPosition(tri_ind,adj_tri_ind,pos_vert +backtraceLength*backtraceDirection,segsIntersection,opposite_edge_ind);
            
            return backtraceFromEdge_woFrame(
                                             dt, opposite_edge_ind, adj_tri_ind, segsIntersection, next_projected_prev_position, previous_positions, weights, triangle_indices, 1 );
        }
        
    }
    
    // When any of the cases 1~3 does not apply. Something went wrong.
    assert( false );
    
    return -1;
}

void HGF::setLocalFrame( const Eigen::Vector3d& localX, size_t tri_ind, Eigen::Matrix3d& end_local_frame )
{
    const Eigen::Vector3d ex = ( localX ).normalized();
    const Eigen::Vector3d ey = triangle_normals[ tri_ind ];
    const Eigen::Vector3d ez = ex.cross( ey );
    
    end_local_frame <<
    ex.x(), ex.y(), ex.z(),
    ey.x(), ey.y(), ey.z(),
    ez.x(), ez.y(), ez.z();
}

int HGF::backtraceFromEdge_withFrame(const double& dt, const size_t e_ind, const size_t tri_ind, const Eigen::Vector3d& intersection, const Eigen::Vector3d& projected_prev_position, std::vector<Eigen::Vector3d>& previous_positions, std::vector<double>& weights, std::vector<size_t>& triangle_indices, std::vector<Eigen::Matrix3d>& end_local_frames, int numCrossedEdges) {
    
#ifdef advect_experiment
    
    if (!intersectionOnce) {
        trace.push_back( intersection );
        // Print out the orbit until now.
        std::cout<<"path lengths\n";
        for(auto path:path_lengths){
            std::cout<<path<<",";
        }std::cout<<"\n\n";
        std::cout<<"trace\n";
        for(auto point:trace){
            std::cout<<point.x()<<","<<point.y()<<","<<point.z()<<"\n";
        }
        std::cout<<"Triangle:\n";
        
        printTri(tri_ind);
    }
#endif //advect_experiment
    using namespace Eigen;
    Vector3d tri_normal = triangle_normals[ tri_ind ];
        
    // In case prev_position is inside the triangle tri_ind
    if (inside_triangle3d( tri_ind, projected_prev_position, 1e-6 )) {//1e-6 standard value
#ifdef advect_experiment
        if (!intersectionOnce) {
            trace.push_back( projected_prev_position );
            path_lengths.push_back( ( projected_prev_position - intersection ).norm() );
            std::cout << "path lengths\n";
            for (auto path : path_lengths) {
                std::cout << path << ",";
            }std::cout << "\n\n";
            std::cout << "trace\n";
            for (auto point : trace) {
                std::cout << point.x() << "," << point.y() << "," << point.z() << "\n";
            }std::cout << "\n";
            
            intersectionOnce = true;
        }
#endif //advect_experiment
        
        // Get the local frame of the triangle
        Eigen::Matrix3d end_local_frame;
        setLocalFrame( projected_prev_position - intersection, tri_ind, end_local_frame );
        
        previous_positions.push_back(projected_prev_position);
        weights.push_back(1);
        triangle_indices.push_back(tri_ind);
        end_local_frames.push_back(end_local_frame);
        
        return 0;
    }
    else if(numCrossedEdges>maxCrossedEdgesInAdvection){
#ifdef advect_experiment
        if (!intersectionOnce) {
            trace.push_back( intersection );
            path_lengths.push_back( 0.0 );
            std::cout << "path lengths\n";
            for (auto path : path_lengths) {
                std::cout << path << ",";
            }std::cout << "\n\n";
            std::cout << "trace\n";
            for (auto point : trace) {
                std::cout << point.x() << "," << point.y() << "," << point.z() << "\n";
            }std::cout << "\n";
            
            intersectionOnce = true;
        }
#endif //advect_experiment
        
        // Get the local frame of the triangle
        Eigen::Matrix3d end_local_frame;
        setLocalFrame( projected_prev_position - intersection, tri_ind, end_local_frame );
        
        previous_positions.push_back(intersection);
        weights.push_back(1);
        triangle_indices.push_back(tri_ind);
        end_local_frames.push_back(end_local_frame);
        
        return 0;
    }
    // In case the trace intersects an edge of the triangle
    else {
        const Vector3d& seg1_begin = intersection;
        const Vector3d& seg1_end = projected_prev_position;
        
        const auto& inc_edges = mesh().m_triangle_to_edge_map[ tri_ind ];
        
        std::unordered_map<size_t, Vector3d> mapEdgeAndIntersections;
        
        for (int ei_t = 0; ei_t < 3; ei_t++) {
            size_t edge = inc_edges[ ei_t ];
            if (edge == e_ind) { continue; }
            
            size_t vi0 = this->mesh().m_edges[ edge ][ 0 ];
            size_t vi1 = this->mesh().m_edges[ edge ][ 1 ];
            
            const auto& v0 = surfTrack()->get_position( vi0 );
            const auto& v1 = surfTrack()->get_position( vi1 );
            
            Vector3d seg2_begin = vc( v0 );
            Vector3d seg2_end = vc( v1 );
            
            Vector3d segsIntersection( 0.0, 0.0, 0.0 );
            
            // See if the orbit intersects an edge of the triangle
            bool isIntersect = twoSegmentsIntersect( segsIntersection, seg1_begin, seg1_end, seg2_begin, seg2_end ) > 0;
            
            if (isIntersect) {
                
                if(surfTrack()->edge_is_all_solid(edge)){
#ifdef advect_experiment
                    if (!intersectionOnce) {
                        trace.push_back(projected_prev_position);
                        path_lengths.push_back((projected_prev_position - intersection).norm());
                        std::cout << "path lengths\n";
                        for (auto path : path_lengths) {
                            std::cout << path << ",";
                        }
                        std::cout << "\n\n";
                        std::cout << "trace\n";
                        for (auto point : trace) {
                            std::cout << point.x() << "," << point.y() << "," << point.z() << "\n";
                        }
                        std::cout << "\n";
                        intersectionOnce = true;
                    }
#endif //advect_experiment
                    
                    Eigen::Matrix3d end_local_frame;
                    // Get the local frame of the triangle
                    setLocalFrame( segsIntersection - intersection, tri_ind, end_local_frame );
                    
                    previous_positions.push_back(segsIntersection);
                    weights.push_back(1);
                    triangle_indices.push_back(tri_ind);
                    end_local_frames.push_back(end_local_frame);
                    
                    return 0;
                }
                mapEdgeAndIntersections[edge]=segsIntersection;
                
            }
        }
        
        if(mapEdgeAndIntersections.size()==0)
        {return -1;}
        
        // Choose edge whose segIntersection is farest. This logic is for avoiding circulation occuring at a point.
        auto comparator=[&intersection]
        (const std::pair<size_t,Vector3d>& edgeAndSegIntersect0,const std::pair<size_t,Vector3d>&edgeAndSegIntersect1)
        {
            return (edgeAndSegIntersect0.second-intersection).squaredNorm()<(edgeAndSegIntersect1.second-intersection).squaredNorm();
            
        };
        const auto& edgeAndSegIntersect=std::max_element(mapEdgeAndIntersections.begin(), mapEdgeAndIntersections.end(),comparator);
        // choose edge with max dist(mapEdgeAndIntersections.begin(), intersection)
        
        const size_t edge=edgeAndSegIntersect->first;
        const Vector3d& segsIntersection=edgeAndSegIntersect->second;
        
#ifdef advect_experiment
        path_lengths.push_back( ( segsIntersection - intersection ).norm() );
#endif //advect_experiment
        
        int backtrace_success = -1;
        
        const bool is_on_Tjunction=mesh().is_edge_nonmanifold(edge);
        std::vector<size_t> valid_adj_triangles;
        std::vector<double> incoming_water_rates;
        
        if (is_on_Tjunction){
            // Finish backtrace if there are no triangles from which water is incoming.
            if(valid_adj_triangles.empty()){
                Eigen::Matrix3d end_local_frame;
                setLocalFrame( projected_prev_position - intersection, tri_ind, end_local_frame );
                end_local_frames.push_back(end_local_frame);
                
                previous_positions.push_back(segsIntersection);
                weights.push_back(1);
                triangle_indices.push_back(tri_ind);
                
#ifdef advect_experiment
                if( !intersectionOnce ) {
                    trace.push_back( segsIntersection );
                    intersectionOnce = true;
                }
#endif //advect_experiment
                return -2;
            }
            
            for(size_t ti_adj=0;ti_adj<valid_adj_triangles.size();++ti_adj){
                const size_t adj_tri_ind=valid_adj_triangles[ti_adj];
                
                // rotate (projected_prev_position - intersection) on the triangle [tri_ind] such that the result is on the triangle [adj_tri_ind].
                Vector3d next_projected_prev_position= findNextProjectedPrevPosition(tri_ind,adj_tri_ind,projected_prev_position,segsIntersection,edge);
                
                std::vector<Vector3d> partial_previous_positions;
                std::vector<double> partial_weights;
                std::vector<size_t> partial_triangle_indices;
                std::vector<Matrix3d> partial_end_local_frames;
                
                backtrace_success = backtraceFromEdge_withFrame(dt, edge, adj_tri_ind, segsIntersection, next_projected_prev_position, partial_previous_positions, partial_weights, partial_triangle_indices, partial_end_local_frames, numCrossedEdges+1);
                
                if(backtrace_success==-1){
                    return -1;
                }
                
                previous_positions.insert(previous_positions.end(), partial_previous_positions.begin(),partial_previous_positions.end());
                triangle_indices.insert(triangle_indices.end(), partial_triangle_indices.begin(),partial_triangle_indices.end());
                end_local_frames.insert(end_local_frames.end(), partial_end_local_frames.begin(),partial_end_local_frames.end());
                
                for(const double& partial_weight: partial_weights){
                    weights.push_back(partial_weight*incoming_water_rates[ti_adj]);
                }
            }
            return -2;
        }else{
            size_t adj_tri_ind = mesh().m_edge_to_triangle_map[ edge ][ 0 ] != tri_ind ? mesh().m_edge_to_triangle_map[ edge ][ 0 ] : mesh().m_edge_to_triangle_map[ edge ][ 1 ];
            // rotate (projected_prev_position - intersection) on the triangle [tri_ind] such that the result is on the triangle [adj_tri_ind].
            
            Vector3d next_projected_prev_position= findNextProjectedPrevPosition(tri_ind,adj_tri_ind,projected_prev_position,segsIntersection,edge);
            return backtraceFromEdge_withFrame(dt, edge, adj_tri_ind, segsIntersection, next_projected_prev_position, previous_positions, weights, triangle_indices, end_local_frames, numCrossedEdges+1);
        }
        
        
    }
    
    return -1;
}

// perform backtrace from edge e_ind after velocity has been determined
int HGF::backtraceFromEdge_woFrame(const double& dt, const size_t e_ind, const size_t tri_ind, const Eigen::Vector3d& intersection, const Eigen::Vector3d& projected_prev_position, std::vector<Eigen::Vector3d>& previous_positions, std::vector<double>& weights, std::vector<size_t>& triangle_indices, int numCrossedEdges) {
    
#ifdef advect_experiment
    if (!intersectionOnce) {
        trace.push_back( intersection );
        // Print out the orbit until now.
        std::cout<<"path lengths\n";
        for(auto path:path_lengths){
            std::cout<<path<<",";
        }std::cout<<"\n\n";
        std::cout<<"trace\n";
        for(auto point:trace){
            std::cout<<point.x()<<","<<point.y()<<","<<point.z()<<"\n";
        }std::cout<<"Triangle:\n";
        
        printTri(tri_ind);
    }
#endif //advect_experiment
    using namespace Eigen;
    const Vector3d& tri_normal = triangle_normals[ tri_ind ];
    
    // In case prev_position is inside the triangle tri_ind
    if (inside_triangle3d( tri_ind, projected_prev_position, 1e-6 )) {
        
#ifdef advect_experiment
        if (!intersectionOnce) {
            trace.push_back( projected_prev_position );
            path_lengths.push_back( ( projected_prev_position - intersection ).norm() );
            std::cout << "path lengths\n";
            for (auto path : path_lengths) {
                std::cout << path << ",";
            }std::cout << "\n\n";
            std::cout << "trace\n";
            for (auto point : trace) {
                std::cout << point.x() << "," << point.y() << "," << point.z() << "\n";
            }std::cout << "\n";
            
            intersectionOnce = true;
        }
#endif //advect_experiment
        previous_positions.push_back(projected_prev_position);
        weights.push_back(1);
        triangle_indices.push_back(tri_ind);
        return 0;
    }
    else if(numCrossedEdges>maxCrossedEdgesInAdvection){
        std::cout<<"Num of crossed edges in advect_vertex reached the limit\n";
#ifdef advect_experiment
        if (!intersectionOnce) {
            trace.push_back( intersection );
            path_lengths.push_back(0.0);
            std::cout << "path lengths\n";
            for (auto path : path_lengths) {
                std::cout << path << ",";
            }std::cout << "\n\n";
            std::cout << "trace\n";
            for (auto point : trace) {
                std::cout << point.x() << "," << point.y() << "," << point.z() << "\n";
            }std::cout << "\n";
            
            intersectionOnce = true;
        }
#endif //advect_experiment
        previous_positions.push_back(intersection);
        weights.push_back(1);
        triangle_indices.push_back(tri_ind);
        return 0;
    }
    // In case the trace intersects an edge of the triangle
    else {
        const Vector3d& seg1_begin = intersection;
        const Vector3d& seg1_end = projected_prev_position;
        
        const auto& inc_edges = mesh().m_triangle_to_edge_map[ tri_ind ];
        
        std::unordered_map<size_t, Vector3d> mapEdgeAndIntersections;
        
        for (int ei_t = 0; ei_t < 3; ei_t++) {
            size_t edge = inc_edges[ ei_t ];
            if (edge == e_ind) { continue; }
            
            size_t vi0 = this->mesh().m_edges[ edge ][ 0 ];
            size_t vi1 = this->mesh().m_edges[ edge ][ 1 ];
            
            const auto& v0 = surfTrack()->get_position( vi0 );
            const auto& v1 = surfTrack()->get_position( vi1 );
            
            Vector3d seg2_begin = vc( v0 );
            Vector3d seg2_end = vc( v1 );
            
            Vector3d segsIntersection( 0.0, 0.0, 0.0 );
            // See if the orbit intersects an edge of the triangle
            bool isIntersect = twoSegmentsIntersect( segsIntersection, seg1_begin, seg1_end, seg2_begin, seg2_end ) > 0;
            
            if (isIntersect) {
                mapEdgeAndIntersections[edge]=segsIntersection;
            }
        }
        
        if(mapEdgeAndIntersections.empty()){return -1;}
        
        // Choose edge whose segIntersection is farthest. This logic is for avoiding circulation occurring at a point.
        auto comparator=[&intersection]
        (const std::pair<size_t,Vector3d>& edgeAndSegIntersect0,const std::pair<size_t,Vector3d>&edgeAndSegIntersect1)
        {
            
            return (edgeAndSegIntersect0.second-intersection).squaredNorm()<(edgeAndSegIntersect1.second-intersection).squaredNorm();
            
        };
        const auto& edgeAndSegIntersect=std::max_element(mapEdgeAndIntersections.begin(), mapEdgeAndIntersections.end(),comparator);
        // choose edge with max dist(mapEdgeAndIntersections.begin(), intersection)
        
        const size_t edge=edgeAndSegIntersect->first;
        const Vector3d& segsIntersection=edgeAndSegIntersect->second;
        
        if( surfTrack()->edge_is_all_solid( edge ) ){
#ifdef advect_experiment
            if( !intersectionOnce ) {
                trace.push_back( segsIntersection );
                intersectionOnce = true;
            }
#endif //advect_experiment
            previous_positions.push_back(segsIntersection);
            weights.push_back(1);
            triangle_indices.push_back(tri_ind);
            
            return 0;
            
        }
        
#ifdef advect_experiment
        path_lengths.push_back( ( segsIntersection - intersection ).norm() );
#endif //advect_experiment
        
        int backtrace_success=-1;
        
        const auto& adj_tris=mesh().m_edge_to_triangle_map[ edge ];
        const bool is_on_Tjunction=mesh().is_edge_nonmanifold(edge);
        std::vector<size_t> valid_adj_triangles;
        std::vector<double> incoming_water_rates;
        if(is_on_Tjunction){
            collectValidTrianglesAndIncomingWaters(dt,tri_ind,edge,valid_adj_triangles, incoming_water_rates);
            
            // Finish backtrace if there are no triangles from which water is incoming.
            if(valid_adj_triangles.empty()){
                previous_positions.push_back(segsIntersection);
                weights.push_back(1);
                triangle_indices.push_back(tri_ind);
                return -2;
            }
            
            for(size_t ti_adj=0;ti_adj<valid_adj_triangles.size();++ti_adj){
                const size_t adj_tri_ind=valid_adj_triangles[ti_adj];
                
                Vector3d next_projected_prev_position= findNextProjectedPrevPosition(tri_ind,adj_tri_ind,projected_prev_position,segsIntersection,edge);
                
                std::vector<Vector3d> partial_previous_positions;
                std::vector<double> partial_weights;
                std::vector<size_t> partial_triangle_indices;
                
                backtrace_success = backtraceFromEdge_woFrame(dt, edge, adj_tri_ind, segsIntersection, next_projected_prev_position, partial_previous_positions, partial_weights, partial_triangle_indices,numCrossedEdges+1 );
                
                if(backtrace_success==-1){
                    return -1;
                }
                
                previous_positions.insert(previous_positions.end(), partial_previous_positions.begin(),partial_previous_positions.end());
                triangle_indices.insert(triangle_indices.end(), partial_triangle_indices.begin(),partial_triangle_indices.end());
                for(const double& partial_weight: partial_weights){
                    weights.push_back(partial_weight*incoming_water_rates[ti_adj]);
                }
            }
            return -2;
            
        }else{
            
            int adj_tri_ind = mesh().m_edge_to_triangle_map[ edge ][ 0 ] != tri_ind ? mesh().m_edge_to_triangle_map[ edge ][ 0 ] : mesh().m_edge_to_triangle_map[ edge ][ 1 ];
            Vector3d next_projected_prev_position= findNextProjectedPrevPosition(tri_ind,adj_tri_ind,projected_prev_position,segsIntersection,edge);
            backtrace_success = backtraceFromEdge_woFrame(dt, edge, adj_tri_ind, segsIntersection, next_projected_prev_position, previous_positions, weights, triangle_indices,numCrossedEdges+1 );
            
            return backtrace_success;
            
        }
        
    }
    
    return -1;
}

void HGF::compute_one_ring_vertices(size_t vertex_index, std::vector<size_t>& adjacent_vertices)
{
    // This function assumes verex is neither boundary nor non-manifold.

    using namespace Eigen;
    // Use them for debugging.
    adjacent_vertices.clear();
    
    if(vertex_index==-1){
        std::cout<<"invalid vertex index"<<std::endl;
        exit(0);
    }
    
    if( (mesh().m_is_boundary_vertex[vertex_index] or mesh().is_vertex_nonmanifold(vertex_index))){
        
        std::cout<<"vertex_index is either boundary or non-manifold in get_1ring_vertices."<<std::endl;
        exit(-1);
    }
    
    // Take one tri
    const auto& inc_tris=mesh().m_vertex_to_triangle_map[vertex_index];
    
    const auto& label=mesh().get_triangle_label(inc_tris[0]);
    for(size_t ti_v=0;ti_v<inc_tris.size();++ti_v){
        const auto& label1=mesh().get_triangle_label(inc_tris[ti_v]);
        if(label1!=label){
            std::cout<<"detect irregular labels around vertex "<<vertex_index<<std::endl;
            mapVerticesIrregular[vertex_index]=true;
            return;
        }
        
    }
    
    const size_t& first_tri=inc_tris[0];
    size_t num_inc_tris=inc_tris.size();
    adjacent_vertices.resize(num_inc_tris);
    
    int next_tri_ind=first_tri;
    const size_t num_e=mesh().ne();
    const size_t num_t=mesh().nt();
    
    for(size_t ti_v=0;ti_v<num_inc_tris;++ti_v){
        const auto& tri=mesh().m_tris[next_tri_ind];
        
        for(int vi_t=0;vi_t<3;++vi_t){
            if(tri[vi_t]==vertex_index){
                size_t prev_vertex= tri[(vi_t+2)%3];
                adjacent_vertices[ti_v]=prev_vertex;
                size_t edge_ind=mesh().get_edge_index(vertex_index, prev_vertex);
                if(edge_ind==num_e){
                    std::cout<<"Edge cannot find in get_1ring_vertices."<<std::endl;
                    exit(-1);
                }
                const auto& tris_e=mesh().m_edge_to_triangle_map[edge_ind];
                if( tris_e.size()!=2){
                    std::cout<<"Edge is either boundary or non-manifold in get_1ring_vertices."<<std::endl;
                    exit(-1);
                }
                next_tri_ind= tris_e[0]==next_tri_ind?tris_e[1]:tris_e[0];
                
                if(next_tri_ind<0 or next_tri_ind>=num_t){
                    std::cout<<"detect irregular triangle configuration around vertex "<<vertex_index<<std::endl;
                    mapVerticesIrregular[vertex_index]=true;
                    return;
                }
            }
        }
        
    }
    
    std::set<size_t> unduplicated_adjacent_vertices;
    for(size_t adj_v:adjacent_vertices){
        unduplicated_adjacent_vertices.insert(adj_v);
    }
    if(adjacent_vertices.size()!=unduplicated_adjacent_vertices.size()){
        std::cout<<"detect irregular labels around vertex "<<vertex_index<<std::endl;
        mapVerticesIrregular[vertex_index]=true;
        
        return;
    }
    
    return;
}


int HGF::semiLagrangianFromTriangle( double dt, size_t tri_ind, Eigen::Vector3d& previous_velocity, std::vector<Eigen::Vector3d>& previous_positions, std::vector<double>& weights, std::vector<size_t>& triangle_indices, std::vector<Eigen::Matrix3d>& end_local_frames) {
    
    using namespace Eigen;
    
#ifdef advect_experiment
    
    if (!intersectionOnce) {
        trace.push_back( triangle_centers[tri_ind] );
        // Print out the orbit until now.
        std::cout<<"path lengths\n";
        for(auto path:path_lengths){
            std::cout<<path<<",";
        }std::cout<<"\n\n";
        std::cout<<"trace\n";
        for(auto point:trace){
            std::cout<<point.x()<<","<<point.y()<<","<<point.z()<<"\n";
        }
        
        std::cout<<"Triangle:\n";
        printTri(tri_ind);
    }
#endif //advect_experiment
    
    Eigen::Vector3d prev_pos=triangle_centers[tri_ind]-dt*previous_velocity;
    
    // If prev_pos is inside triangle, return as is.
    // If orbit intersects an edge, then return trackTriangle_edge.
    if (inside_triangle3d( tri_ind, prev_pos, 1e-6 )) {//1e-6 standard value
#ifdef advect_experiment
        if (!intersectionOnce) {
            trace.push_back( prev_pos );
            path_lengths.push_back( ( prev_pos - triangle_centers[tri_ind] ).norm() );
            std::cout << "path lengths\n";
            for (auto path : path_lengths) {
                std::cout << path << ",";
            }std::cout << "\n\n";
            std::cout << "trace\n";
            for (auto point : trace) {
                std::cout << point.x() << "," << point.y() << "," << point.z() << "\n";
            }std::cout << "\n";
            
            intersectionOnce = true;
        }
#endif //advect_experiment
        
        // Get local frame of the triangle.
        // In this case end_local_frame == begin_local_frame
        // end_local_frame is the return value
        // push back new end_local_frame i.e. Matrix3d to end_local_frames
        // and do the below for end_local_frames[0]
        // set triangle_indices, previous_positions, and weights
        Eigen::Matrix3d end_local_frame;
        setLocalFrame( prev_pos - triangle_centers[ tri_ind ] , tri_ind, end_local_frame );
        
        previous_positions.push_back(prev_pos);
        weights.push_back(1);
        triangle_indices.push_back(tri_ind);
        end_local_frames.push_back(end_local_frame);
        
        return 0;
        
        // TODO: store previous position(s), weights(1), triangle_index(s), local_frame(s)
    }
    else
    {
        
        const Vector3d& seg1_begin = triangle_centers[tri_ind];
        const Vector3d& seg1_end = prev_pos;
        const auto& inc_edges = mesh().m_triangle_to_edge_map[ tri_ind ];
        
        for (int ei_t = 0; ei_t < 3; ei_t++) {
            size_t edge = inc_edges[ ei_t ];
            
            size_t vi0 = this->mesh().m_edges[ edge ][ 0 ];
            size_t vi1 = this->mesh().m_edges[ edge ][ 1 ];
            
            const auto& v0 = surfTrack()->get_position( vi0 );
            const auto& v1 = surfTrack()->get_position( vi1 );
            
            Vector3d seg2_begin = vc( v0 );
            Vector3d seg2_end = vc( v1 );
            
            Vector3d segsIntersection( 0.0, 0.0, 0.0 );
            bool isIntersect = twoSegmentsIntersect( segsIntersection, seg1_begin, seg1_end, seg2_begin, seg2_end ) > 0;
            
            if (isIntersect) {
#ifdef advect_experiment
                path_lengths.push_back( ( segsIntersection - triangle_centers[tri_ind] ).norm() );
#endif //advect_experiment
                
                // if edge is boundary then, treat this inside this triangle.
                if(surfTrack()->edge_is_all_solid(edge)){
#ifdef advect_experiment
                    if( !intersectionOnce ) {
                        trace.push_back( segsIntersection );
                        intersectionOnce = true;
                    }
#endif //advect_experiment
                    
                    prev_pos=segsIntersection;
                    
                    //  store all the information needed i.e. prev_positions,weights, tri_indices, local_frame
                    // finally return 0 instead of tri_ind
                    
                    // Get local frame of the triangle.
                    // In this case end_local_frame == begin_local_frame
                    Eigen::Matrix3d end_local_frame;
                    setLocalFrame( prev_pos - triangle_centers[ tri_ind ], tri_ind, end_local_frame );
                    
                    previous_positions.push_back(prev_pos);
                    weights.push_back(1);
                    triangle_indices.push_back(tri_ind);
                    end_local_frames.push_back(end_local_frame);
                    
                    return 0;
                    
                }
                
                int backtrace_success = -1;
                // If we are on a T-junction, decide if we should split path for each adjacent triangles.
                
                // We can use collectValidTrianglesAndIncomingWaters functions.
                // just run through valid_adj_triangles
                const bool is_on_Tjunction=mesh().is_edge_nonmanifold(edge);
                std::vector<size_t> valid_adj_triangles;
                std::vector<double> incoming_water_rates;
                
                if(is_on_Tjunction) {
#ifdef advect_experiment
                    if (!intersectionOnce) {
                        trace.push_back(segsIntersection);
                        intersectionOnce = true;
                    }
#endif //advect_experiment
                    collectValidTrianglesAndIncomingWaters(dt,tri_ind,edge,valid_adj_triangles, incoming_water_rates);
                    
                    // Finish backtrace if there are no triangles from which water is incoming.
                    if (valid_adj_triangles.empty()) {
                        
                        previous_positions.push_back(segsIntersection);
                        weights.push_back(1);
                        triangle_indices.push_back(tri_ind);
                        Eigen::Matrix3d end_local_frame;
                        setLocalFrame(prev_pos - triangle_centers[ tri_ind ], tri_ind, end_local_frame);
                        end_local_frames.push_back(end_local_frame);
                        
                        return -2;
                    }
                    
                    for(size_t ti_adj=0;ti_adj<valid_adj_triangles.size();++ti_adj){
                        const size_t adj_tri_ind=valid_adj_triangles[ti_adj];
                        
                        Vector3d next_projected_prev_position=
                        findNextProjectedPrevPosition(tri_ind,adj_tri_ind,prev_pos,segsIntersection,edge);
                        
                        std::vector<Vector3d> partial_previous_positions;
                        std::vector<double> partial_weights;
                        std::vector<size_t> partial_triangle_indices;
                        std::vector<Matrix3d> partial_end_local_frames;
                        
                        backtrace_success = backtraceFromEdge_withFrame(dt, edge, adj_tri_ind, segsIntersection, next_projected_prev_position, partial_previous_positions, partial_weights, partial_triangle_indices, partial_end_local_frames, 1);
                        
                        if(backtrace_success == -1){
                            return -1;
                        }
                        
                        previous_positions.insert(previous_positions.end(), partial_previous_positions.begin(),partial_previous_positions.end());
                        triangle_indices.insert(triangle_indices.end(), partial_triangle_indices.begin(),partial_triangle_indices.end());
                        end_local_frames.insert(end_local_frames.end(), partial_end_local_frames.begin(),partial_end_local_frames.end());
                        for(const double& partial_weight: partial_weights){
                            weights.push_back(partial_weight*incoming_water_rates[ti_adj]);
                        }
                        
                    }
                    return -2;
                }else{
                    const auto& adj_tris=mesh().m_edge_to_triangle_map[ edge ];
                    const size_t adj_tri_ind=adj_tris[0]!=tri_ind?adj_tris[0]:adj_tris[1];
                    
                    Vector3d next_projected_prev_position=
                    findNextProjectedPrevPosition(tri_ind,adj_tri_ind,prev_pos,segsIntersection,edge);
                    
                    return backtraceFromEdge_withFrame(dt, edge, adj_tri_ind, segsIntersection, next_projected_prev_position, previous_positions, weights, triangle_indices, end_local_frames, 1);
                }
                
            }
        }
        
    }
    
    return -1;
}

// this fuction performs the backtrace; contains function calcVertexVelocity to determine the backtrace velocity
// input: time step, source index, output vectors of backtraced positions, weights (due to path splitting), trinagle indices
//     (for local coordinates of backtraced positions), magnitude of backtrace velocity, iteration of backtrace
int HGF::semiLagrangianFromVertex(double dt, size_t v_ind,std::vector<Eigen::Vector3d>& previous_positions, std::vector<double>& weights, std::vector<size_t>& triangle_indices, double magPrevVelocity, int ite){
    using namespace Eigen;
    
    // list of triangles adjacent to source vertex
    const auto& inc_tris = mesh().m_vertex_to_triangle_map[v_ind];
    const size_t num_inc_tris = inc_tris.size();
    
    // average of velocities on tris adjacent to v_ind when flattened onto plane
    Vector3d ave_tri_velocities2d;
    std::unordered_map <size_t,Eigen::Vector3d> modified_ad_positions;
    
    // Calculate the backtrace velocity
    calcVertexVelocity(v_ind, -1, ave_tri_velocities2d, MatrixXd(),modified_ad_positions);
    
    // Randomly pertubate direction later than the first iteration to avoid failure of backtrace.
    if(ite>0) {
        const double minorModification=std::pow(0.1, (ite%6+3) );
        ave_tri_velocities2d[ite%3]*=(1.0-ite*minorModification);
    }
    
    // in case the projected_prev_velocity2d is nearly 0, set the source vertex position as the backtraced position
    // and save first incident triangle; also return this triangle's index
    const double eps = 1e-8;
    if( ave_tri_velocities2d.norm() < eps ) {
        previous_positions.push_back(pos( v_ind ));
        weights.push_back(1);
        const size_t triangle_index = mesh().m_vertex_to_triangle_map[ v_ind ][ 0 ];
        triangle_indices.push_back(triangle_index);
        return triangle_index;
    }
    
    // Set backtrace velocity to magPrevVelocity.
    // Can only evaluate as false if not given (magPrevVelocity is an optional parameter)
    // magPrevVelocity only passed when vertex velocity is taken from LosTopos
    if(magPrevVelocity>-1.0) {
        ave_tri_velocities2d.normalize();
        ave_tri_velocities2d*=magPrevVelocity;
    }
    
    // Naive previous position in local coordinates.
    // Essentially just a scaling of the velocity on a vertex by the time step size
    Vector3d prev_position2d = -dt * ave_tri_velocities2d;
    
    // Iterate over the incident triangles
    for( size_t ti = 0; ti < num_inc_tris; ++ti ) {
        
        // For debug purpose
#ifdef advect_experiment
        if( !intersectionOnce ) {
            trace.clear();
            path_lengths.clear();
        }
#endif //advect_experiment
        
        size_t src_tri = inc_tris[ ti ];
        
        // Save naive previous position.
        // Temporary set previous position to be prev_position2d.
        
        // Backtrace is done inside this function
        // Return -2 if it hits T-junction edge
        
        // Success is 0 if it found a single final destination, -1 if backtrace failed, -2 if it hits t-junction edge (can be successfull too).
        int success = backtraceFromVertex( dt, v_ind, src_tri, ave_tri_velocities2d, modified_ad_positions, prev_position2d, previous_positions, weights, triangle_indices );
        
        if( success != -1 ) {
            
            return success;
        }
    }
    
    return -1;
    
}

double HGF::advect_triangle_BFECC(double dt, const Eigen::MatrixXd& AdvectionSource, Eigen::MatrixXd& AdvectionTarget)
{
    using namespace Eigen;
    
    bool local_energy_preserve=true;
    MatrixXd vertexVelocity=AdvectionSource;
    MatrixXd AdvectionTarget0=AdvectionSource;
    
    advect_triangle_standard(dt, vertexVelocity, AdvectionTarget0);
    // Update vertexVelocity for next advection
    vectorTri2vert(AdvectionTarget0, vertexVelocity, local_energy_preserve);
    
    advect_triangle_standard(-dt, vertexVelocity, AdvectionTarget0);
    AdvectionTarget.noalias()=1.5*AdvectionTarget-0.5*AdvectionTarget0;
    // Update vertexVelocity for next advection
    vectorTri2vert(AdvectionTarget, vertexVelocity, local_energy_preserve);
    
    advect_triangle_standard(dt, vertexVelocity, AdvectionTarget);
    
    return dt;
}

double HGF::advect_triangle_MacCormack(double dt, const Eigen::MatrixXd& AdvectionSource, Eigen::MatrixXd& AdvectionTarget)
{
    using namespace Eigen;
    
    // not shi and yu version
    bool local_energy_preserve=true;
    MatrixXd vertexVelocity=AdvectionSource;
    MatrixXd AdvectionTarget0=AdvectionTarget;
    
    advect_triangle_standard(dt, vertexVelocity, AdvectionTarget0);
    // Update vertexVelocity for next advection
    vectorTri2vert(AdvectionTarget0, vertexVelocity, local_energy_preserve);
    // Save first advection data
    MatrixXd AdvectionTargetAfterFirstAdvect=AdvectionTarget0;
    
    advect_triangle_standard(-dt, vertexVelocity, AdvectionTarget0);
    AdvectionTarget.noalias()=AdvectionTargetAfterFirstAdvect + 0.5*(AdvectionTarget-AdvectionTarget0);
    
    
    return dt;
}

double HGF::advect_triangleVelocity(double dt, TriangleAdvectType advectType){
    using namespace Eigen;
    
    // FIXME: Do not use this if we are in simulation mode 2 (flow+deformation), and proceed vertex velocity -> triangle velocity -> advection -> projection -> vertex velocity
    // since we already have vertex velocity in the beginning of the time step.
    velocityTriangle2Vertex(true);
    
    const MatrixXd sourceMat=fvmat();
    MatrixXd triVelMat=getTriVelocity();
    
    switch(advectType){
        case TriStandard:
            advect_triangle_standard(dt, sourceMat, triVelMat);
            break;
        case TriBFECC:
            advect_triangle_BFECC(dt, sourceMat, triVelMat);
            break;
        case TriMacCormack:
            advect_triangle_MacCormack(dt, sourceMat, triVelMat);
            break;
    }
    
    setTriVelocity(triVelMat);
    
    return dt;
}

// In most cases, AdvectionSource is num_t x 3 matrices that store vertex velocities.
// AdvectionTarget is the output, which is again num_t x matrices storing triangle velocities.
double HGF::advect_triangle_standard(double dt, const Eigen::MatrixXd& AdvectionSource, Eigen::MatrixXd& AdvectionTarget) {
    
    using namespace Eigen;
    const bool forScalar = AdvectionSource.cols() == 1;
    
    size_t num_t = mesh().nt();
    
    AdvectionTarget.resize(num_t, AdvectionSource.cols());
    
    std::vector<Vector3d> new_tri_velocities = triangle_velocities;
    std::vector<double> new_tri_values(num_t);
    
    std::unordered_map<size_t,bool> trianglesDestinationUnfound;
    
    //for( size_t ti = 0; ti < num_t; ++ti ){
    for (size_t ti : real_triangles) {
        Vector3d velocity = Vector3d::Zero();
        double newValue = 0.0;
        
        // Local frame of the triangle where the starting point lies.
        Eigen::Matrix3d begin_local_frame = Eigen::Matrix3d::Zero();
        // Local frame of the triangle where the orbit arrives.
        Eigen::Matrix3d end_local_frame = Eigen::Matrix3d::Zero();
        
        // 0: success, -1: failure, -2: success but crossed a T-junction edge.
        int backtrace_success = -1;
        int ite = 0, max_iteration = 30;
        
        Eigen::Vector3d previous_velocity;
        std::vector<Eigen::Vector3d> previous_positions;
        std::vector<double> weights;
        std::vector<size_t> triangle_indices;
        std::vector<Eigen::Matrix3d> end_local_frames;
        
        while (ite < max_iteration and backtrace_success == -1) {
            previous_velocity = triangle_velocities[ ti ];
            if (ite > 0) {
                const double minorModification = std::pow(0.1, ((4 - ite) % 4 + 1));
                previous_velocity[ite % 3] *= (1.0 - ite * minorModification);
                //prev_velocity = Quaterniond( AngleAxisd( minorModification* ite , Vector3d( ite % 3, ite % 3 + 1, ite % 3 + 2 ) ) ) *prev_velocity;
                projectToTangent(previous_velocity, triangle_normals[ti]);
                printVec(previous_velocity);
            }
            
            // Local frame of the triangle where the starting point lies.
            Eigen::Vector3d prev_pos = triangle_centers[ ti ] - dt * previous_velocity;
            Eigen::Vector3d bx = (prev_pos - triangle_centers[ ti ]).normalized();
            Eigen::Vector3d by = (triangle_normals[ ti ]).normalized();
            Eigen::Vector3d bz = (bx.cross(by)).normalized();
            
            begin_local_frame <<
            bx.x(), bx.y(), bx.z(),
            by.x(), by.y(), by.z(),
            bz.x(), bz.y(), bz.z();
            
            backtrace_success = semiLagrangianFromTriangle(dt, ti, previous_velocity, previous_positions, weights, triangle_indices, end_local_frames);
            
            if (backtrace_success == -1 and ite % 10 == 0) {
                std::cout << ite << " times finding tri_ind for triangle " << ti << ".\n";
            }
            ++ite;
        }
        
        if (backtrace_success == -1) {
            std::cerr << "Could not find target triangle for triangle "<<ti <<" in advect_triangleVelocity after " << max_iteration <<
            " times of backtraces." << std::endl;
            
            trianglesDestinationUnfound[ti]=true;
            continue;
        }
        
        if (forScalar) {
            // Currently it is disabled as we are not advecting scalar values on triangles.
            // If you want to use this, adjust the next line properly to the current implementation.
            //interpolateScalarFromVertices(tri_ind, naive_prev_position, AdvectionSource, newValue);
            new_tri_values[ ti ] = newValue;
        } else {
            
            bool local_energy_preserve = 0;
            
            //  for the time being, just use this with tri_indices[0],prev_positions[0],.. and test it with inputmesh_sphere.txt
            // triangle_indices, previous_positions, weights come from backtracing, AdvectionSource is num_t x 3 matrix that stores vertex velocities
            // velocity is the output value, local_energy_preserve flags an energy preservation method of interpolation
            interpolateVectorFromVerticesOnWeightedPaths(triangle_indices, previous_positions, weights, begin_local_frame, end_local_frames, AdvectionSource, velocity, local_energy_preserve);
            
            const bool length_preserve = true;
            
            projectToTangent(velocity, triangle_normals[ ti ], length_preserve);
            new_tri_velocities[ ti ] = velocity;
            
        }
        
    }
    
    // Interpolate velocities of triangles that have not found destinations from their neighbor triangles.
    for(auto triAndUnfound:trianglesDestinationUnfound){
        
        size_t tri_wo_dest=triAndUnfound.first;
        
        std::vector<size_t> adj_tris;
        mesh().get_adjacent_triangles(tri_wo_dest, adj_tris);
        
        size_t num_valid_adj_tris=0;
        new_tri_velocities[tri_wo_dest].setZero();
        
        for(size_t adj_tri:adj_tris){
            if(not trianglesDestinationUnfound[adj_tri] and triangles_to_real[adj_tri]!=-1){
                new_tri_velocities[tri_wo_dest]+=new_tri_velocities[adj_tri];
                ++num_valid_adj_tris;
            }
            if(num_valid_adj_tris>0){
                new_tri_velocities[tri_wo_dest]/=num_valid_adj_tris;
            }
            
        }
        
    }
    
    if (forScalar) {
        for (size_t ti = 0; ti < num_t; ++ti)
            AdvectionTarget(ti, 0) = new_tri_values[ ti ];
    } else {
        for (size_t ti = 0; ti < num_t; ++ti) {
            AdvectionTarget.row(ti) = new_tri_velocities[ ti ];
        }
        
    }
    
    return dt;
}

double HGF::advect_triangle_viaFlux(double dt,Eigen::VectorXd& AdvectionTarget, const Eigen::MatrixXd& TriVelMat){
    
    using namespace Eigen;
    
    const MatrixXd& TriVelocity=TriVelMat.size()==0?getTriVelocity():TriVelMat;
    
    const size_t num_t=mesh().nt();
    const size_t num_e=mesh().ne();
    const size_t num_v=mesh().nv();
    
    VectorXd PrevTriangleScalar=AdvectionTarget;
    
    // Perform mass-preserving advection
    for(size_t ei=0;ei<num_e;++ei){
        const auto tris=mesh().m_edge_to_triangle_map[ei];

        const size_t num_inc_tris=tris.size();
        
        for(size_t t0=0;t0<num_inc_tris;++t0){
            const size_t tri_ind0=tris[t0];
            if(triangles_to_real[tri_ind0]==-1){
                continue;
            }
            for(size_t t1=t0+1;t1<num_inc_tris;++t1){
                const size_t tri_ind1=tris[t1];
                if(triangles_to_real[tri_ind1]==-1){
                    continue;
                }
                
                const size_t va=mesh().m_edges[ei][0];
                const size_t vb=mesh().m_edges[ei][1];
                
                const Vector3d edge = pos(vb)-pos(va );
                
                Vector3d edgeNormalUnNormalized0to1=(triangle_centers[tri_ind1]-triangle_centers[tri_ind0]).normalized()*edge.norm();
                if(edgeNormalUnNormalized0to1.dot(edge_centers[ei]-triangle_centers[tri_ind0])<0){
                    edgeNormalUnNormalized0to1*=-1;
                }
                
                Vector3d edgeVelocity =0.5*(triangle_areas[tri_ind0]* TriVelocity.row(tri_ind0)+triangle_areas[tri_ind1]* TriVelocity.row(tri_ind1))/(triangle_areas[tri_ind0]+triangle_areas[tri_ind1]);
                // Projected to edge normal.
                //edgeVelocity=edgeVelocity.dot(edgeNormalUnNormalized0to1)*edgeNormalUnNormalized0to1();
                
                double value=edgeVelocity.dot(edgeNormalUnNormalized0to1)>=0?PrevTriangleScalar[tri_ind0]:PrevTriangleScalar[tri_ind1];
                
                double delta = value*dt*edgeVelocity.dot(edgeNormalUnNormalized0to1);
                AdvectionTarget[tri_ind0]-=delta/triangle_areas[tri_ind0];
                AdvectionTarget[tri_ind1]+=delta/triangle_areas[tri_ind1];
            }
        }
        
    }
    return dt;
}

double HGF::advect_vertex(double dt, Eigen::VectorXd& AdvectionTarget, VertexAdvectType advectType){
    switch (advectType) {
        case VertStandard:
            advect_vertex_standard(dt, AdvectionTarget);
            break;
        case VertBFECC:
            advect_vertex_BFECC(dt, AdvectionTarget);
            break;
        case VertMacCormack:
            advect_vertex_MacCormack(dt, AdvectionTarget);
            break;
            
    }
    return dt;
}

// advect values on vertices (thickness, derivatives of thickness, color,...)
// returns just value of dt, AdvectionTarget is globally defined and contains the advected values
double HGF::advect_vertex_standard(double dt, Eigen::VectorXd& AdvectionTarget){
    
    std::cout<<"Advect vertex starts\n";
        
    // First perform advection for vertices not adjacent to non-manifold junctions while detecting vertices near non-manifold junctions.
    std::vector<size_t> vertices_near_tjunctions;
    std::vector<size_t> anomalous_vertices;
    
    advect_vertex_core(dt, AdvectionTarget,  &vertices_near_tjunctions,&anomalous_vertices);
    
    const size_t num_v=mesh().nv();
    // For vertices near T-junctions, assign values computed by flux-based approach.
    std::vector<size_t> t_junctions;
    for(size_t vi=0;vi<num_v;++vi){
        if(mesh().is_vertex_nonmanifold(vi)){
            t_junctions.push_back(vi);
        }
    }
    t_junctions.insert(t_junctions.end(),anomalous_vertices.begin(),anomalous_vertices.end());
    
    vertices_near_tjunctions.insert(vertices_near_tjunctions.end(), anomalous_vertices.begin(),anomalous_vertices.end());
    
    const std::vector<size_t>& vertices_for_flux_advection=t_junctions;
    
    // Peform flux-based advection for vertices on junctions.
    if(not vertices_for_flux_advection.empty()){
        Eigen::VectorXd FluxBasedAdvectionTarget=AdvectionTarget;
        advect_vertex_viaTriangle(dt, FluxBasedAdvectionTarget);
        
        for(const auto& vert:vertices_for_flux_advection) {
            AdvectionTarget[vert]=FluxBasedAdvectionTarget[vert];
        }
        
        // For vertices near junctions, blend advection results with neighbor vertices.
        const double smoothingCoef=0.5;
        if(smoothingCoef!=0.0){
            const size_t num_v=mesh().nv();
            
            enum SmoothingTarget{
                Tjunctions,
                VerticesNearTjunctions
            }smoothingTarget=VerticesNearTjunctions;
            
            std::vector<size_t> vertices_for_smoothing;
            switch(smoothingTarget){
                case Tjunctions:{
                    vertices_for_smoothing=t_junctions;
                    break;
                }
                case VerticesNearTjunctions:{
                    vertices_for_smoothing=vertices_near_tjunctions;
                    break;
                }
            }
            
            for(size_t vi_s=0;vi_s<vertices_for_smoothing.size();++vi_s){
                const size_t vert=vertices_for_smoothing[vi_s];
                
                std::vector<size_t> adj_verts;
                mesh().get_adjacent_vertices(vert, adj_verts);
                
                double neighborSum=0.0;
                size_t num_manifold_neighbors=0;
                for(const size_t& adj_v:adj_verts){
                    if(not mesh().is_vertex_nonmanifold(adj_v)){
                        neighborSum+=FluxBasedAdvectionTarget[adj_v];
                        ++num_manifold_neighbors;
                    }
                }
                
                if(num_manifold_neighbors>0){ AdvectionTarget[vert]=(1.0-smoothingCoef)*FluxBasedAdvectionTarget[vert]+smoothingCoef*neighborSum/num_manifold_neighbors;
                }
            }
        }
    }
    
    std::cout<<"Advection on vertices finished.\n";
    return dt;
}

double HGF::advect_vertex_viaTriangle(double dt, Eigen::VectorXd& AdvectionTarget){

    using namespace Eigen;
    
    // Define scalar values on triangles from vertices.
    VectorXd TriangleScalar;
    scalarVert2triNaive(AdvectionTarget, TriangleScalar);
    
    // Peform flux-based advection for triangles
    VectorXd PrevTriangleScalar=TriangleScalar;
    advect_triangle_viaFlux(dt, TriangleScalar);
    
    // Retrieve values on vertices.
    scalarTri2vertNaive(TriangleScalar,AdvectionTarget);
    
    return dt;
}

//advect values on vertices (thickness)
double HGF::advect_vertex_core( double dt, Eigen::VectorXd& AdvectionTarget, std::vector<size_t>* vertices_near_tjunctions,std::vector<size_t>* anomalous_vertices )
{
    
    using namespace Eigen;
    size_t num_v = mesh().nv();
    
    std::vector<Vector3d> new_vert_velocities( num_v );
    std::vector<double> new_vert_values(AdvectionTarget.data(), AdvectionTarget.data()+num_v);
    
    const MatrixXd fvMat = fvmat();
    
    for( size_t vi : non_solid_vertices ) {
        
        if(mesh().is_vertex_nonmanifold(vi)){
            
            // add vi to vertices_near_Tjunctions
            if(vertices_near_tjunctions)
            {vertices_near_tjunctions->push_back(vi);}
            continue;
        }
        
        if(mapVerticesIrregular[vi] and anomalous_vertices){
            std::cout<<"\nanomaly:"<<vi<<" in advection.\n";
            anomalous_vertices->push_back(vi);
            const auto& inc_tris=mesh().m_vertex_to_triangle_map[vi];
            
            continue;
        }
        
        int backtrace_success = -1;
        int ite=0,max_iteration = 30;
        Vector3d prev_position;
        
        std::vector<Eigen::Vector3d> previous_positions;
        std::vector<double> weights;
        std::vector<size_t> triangle_indices;
        
        while(ite<max_iteration and backtrace_success==-1)
        {
            // Backtrace is done here.
            // Here we want to know if the vertex is near T-junctions i.e. if the path hits a T-junction edge or not.
            
            // Perform semi-Lagrangian backtrace from a vertex for given: time step, source vertex index, vector of backtraced positions, weights and
            //  triangle indices, and current iteration.
            
            // return values: -1 = backtrace failed, -2 = we hit a non-mfld edge
            backtrace_success = semiLagrangianFromVertex( dt, vi, previous_positions, weights, triangle_indices, -1.0,ite);
            
            if (vertices_near_tjunctions and backtrace_success ==-2) {
                vertices_near_tjunctions->push_back(vi);
            }
            
            if (backtrace_success==-1 and  ite % 10 == 0) {
                std::cout <<"Backtrace " << ite << " times finding tri_ind for vertex " << vi  <<".\n";
                printVec(pos(vi));
            }++ite;
        }
                
        if(backtrace_success==-1 and anomalous_vertices){
            std::cerr<<"Could not find target triangle for vertex "<<vi<<" in advect_vertex after "<<max_iteration<<
            " times of backtraces."<<std::endl;
            anomalous_vertices->push_back(vi);
            
        }
        
        // Perform interpolation based on previous_positions for given: weights, backtraced positions, triangles of the positions,
        //     current values on the vertices and new values on the vertices
        InterpolateScalarFromVerticesOfMultipleTriangles(weights,previous_positions, triangle_indices, AdvectionTarget, new_vert_values[ vi ]);
        
    }
    
    // Copy new values on vertices to the output vector
    for( size_t vi = 0; vi < num_v; ++vi )
        AdvectionTarget[ vi ] = new_vert_values[ vi ];
    
    return dt;
}

double HGF::advect_vertex_BFECC( double dt, Eigen::VectorXd& AdvectionTarget)
{
    using namespace Eigen;
    const bool limiter=true;
    
    const size_t num_v=mesh().nv();
    VectorXd AdvectionTargetInit;
    std::vector<std::vector<size_t>> prev_triangle_indices(num_v);
    std::vector<std::vector<double>> global_weights(num_v);
    
    std::unordered_map<size_t,bool> mapVerticesIrregularForLimiter;
    if(limiter){
        AdvectionTargetInit=AdvectionTarget;
        
        // for irregular_vertices, we do not use limiter.
        for( size_t vi:real_vertices ) {
            if(surfTrack()->vertex_is_any_solid(vi)
               or mesh().is_vertex_nonmanifold(vi)
               or mapVerticesIrregular[vi]
               )
            {
                continue;
            }
            
            int backtrace_success = -1;
            int ite=0,max_iteration = 30;
            while(ite<max_iteration and backtrace_success==-1)
            {
                
                std::vector<Eigen::Vector3d> previous_positions;
                
                backtrace_success=semiLagrangianFromVertex( dt, vi, previous_positions, global_weights[vi], prev_triangle_indices[vi],-1.0,ite);
                
                if (backtrace_success==-1 and  ite % 10 == 0) {
                    std::cout << ite << " times finding tri_ind for vertex" << vi << ".\n";
                }++ite;
            }
            if(backtrace_success==-1){
                std::cerr<<"Could not find target triangle in advect_vertex_BFECC's semiLagrangian after "<<max_iteration<<
                " times of backtraces. Set "<<vi<<"th vertex as an irregular vertex." <<std::endl;
                
                mapVerticesIrregularForLimiter[vi]=true;
            }
            
        }
    }
    
    VectorXd AdvectionTarget0=AdvectionTarget;
    advect_vertex_standard(dt, AdvectionTarget0);
    advect_vertex_standard(-dt, AdvectionTarget0);
    AdvectionTarget.noalias()=1.5*AdvectionTarget-0.5*AdvectionTarget0;
    advect_vertex_standard(dt, AdvectionTarget);
    
    if(limiter){
        for( const size_t& vi : real_vertices ) {
            //for( size_t vi = 0; vi < num_v; ++vi ) {
            
            if(surfTrack()->vertex_is_any_solid(vi) or mesh().is_vertex_nonmanifold(vi) or mapVerticesIrregular[vi] or mapVerticesIrregularForLimiter[vi]){
                continue;
            }
            
            double max_value=0.0;
            double min_value=0.0;
            
            const std::vector<size_t>&prev_tri_inds=prev_triangle_indices[vi];
            
            for(size_t ti_v=0;ti_v<prev_tri_inds.size();++ti_v){
                
                const auto& prevTri=mesh().m_tris[prev_tri_inds[ti_v]];
                
                double local_max_value = AdvectionTargetInit[prevTri[0]];
                double local_min_value = AdvectionTargetInit[prevTri[0]];
                
                for(int i=1; i<3; i++) {
                    if(local_max_value < AdvectionTargetInit[prevTri[i]]) local_max_value = AdvectionTargetInit[prevTri[i]];
                    if(local_min_value > AdvectionTargetInit[prevTri[i]]) local_min_value = AdvectionTargetInit[prevTri[i]];
                }
                
                max_value+=local_max_value*global_weights[vi][ti_v];
                min_value+=local_min_value*global_weights[vi][ti_v];
                
            }
            
            AdvectionTarget[vi] = clamp(AdvectionTarget[vi], min_value, max_value);
        }
    }
    return dt;
}

double HGF::advect_vertex_MacCormack( double dt, Eigen::VectorXd& AdvectionTarget)
{
    using namespace Eigen;
    const bool limiter=true;
    const size_t num_v=mesh().nv();
    VectorXd AdvectionTargetInit;
    std::vector<std::vector<size_t>> prev_triangle_indices(num_v);
    std::vector<std::vector<double>> global_weights(num_v);
    
    std::unordered_map<size_t,bool> mapVerticesIrregularForLimiter;
    
    if(limiter){
        AdvectionTargetInit=AdvectionTarget;
        
        // for irregular_vertices, we do not use limiter.
        for( size_t vi :real_vertices ) {
            if(surfTrack()->vertex_is_any_solid(vi) or mesh().is_vertex_nonmanifold(vi) or mapVerticesIrregular[vi]){
                continue;
            }
            
            int backtrace_success = -1;
            int ite=0,max_iteration = 30;
            while(ite<max_iteration and backtrace_success==-1)
            {
                
                std::vector<Eigen::Vector3d> previous_positions;
                
                backtrace_success=semiLagrangianFromVertex( dt, vi, previous_positions, global_weights[vi], prev_triangle_indices[vi],-1.0,ite);
                
                if (backtrace_success==-1 and  ite % 10 == 0) {
                    std::cout << ite << " times finding tri_ind for vertex" << vi << ".\n";
                }++ite;
            }
            if(backtrace_success==-1){
                
                std::cerr<<"Could not find target triangle in advect_vertex_Maccormark's semiLagrangian after "<<max_iteration<<
                " times of backtraces. Set "<<vi<<"th vertex as an irregular vertex." <<std::endl;
                
                mapVerticesIrregularForLimiter[vi]=true;
            }
            
        }
    }
    
    VectorXd AdvectionTarget0=AdvectionTarget;
    advect_vertex_standard(dt, AdvectionTarget0);
    VectorXd AdvectionTargetAfterFirstAdvect=AdvectionTarget0;
    advect_vertex_standard(-dt, AdvectionTarget0);
    AdvectionTarget.noalias()=AdvectionTargetAfterFirstAdvect + 0.5*(AdvectionTarget-AdvectionTarget0);
    
    if(limiter){
        for( const size_t& vi : real_vertices ) {
            //for( size_t vi = 0; vi < num_v; ++vi ) {
            
            if(surfTrack()->vertex_is_any_solid(vi) or mesh().is_vertex_nonmanifold(vi) or mapVerticesIrregular[vi] or mapVerticesIrregularForLimiter[vi]){
                continue;
            }
            
            double max_value=0.0;
            double min_value=0.0;
            
            const std::vector<size_t>&prev_tri_inds=prev_triangle_indices[vi];
            
            for(size_t ti_v=0;ti_v<prev_tri_inds.size();++ti_v){
                
                const auto& prevTri=mesh().m_tris[prev_tri_inds[ti_v]];
                
                double local_max_value = AdvectionTargetInit[prevTri[0]];
                double local_min_value = AdvectionTargetInit[prevTri[0]];
                
                for(int i=1; i<3; i++) {
                    if(local_max_value < AdvectionTargetInit[prevTri[i]]) local_max_value = AdvectionTargetInit[prevTri[i]];
                    if(local_min_value > AdvectionTargetInit[prevTri[i]]) local_min_value = AdvectionTargetInit[prevTri[i]];
                }
                
                max_value+=local_max_value*global_weights[vi][ti_v];
                min_value+=local_min_value*global_weights[vi][ti_v];
                
            }
            
            AdvectionTarget[vi] = clamp(AdvectionTarget[vi], min_value, max_value);
        }
    }
    
    return dt;
}

void HGF::project_flow_v_tangent() {
    int num_v = mesh().nv();
    
    for (size_t vi = 0; vi < num_v; ++vi) {
        if(mesh().is_vertex_nonmanifold(vi)){continue;}
        auto vertex_normal = vc( surfTrack()->get_vertex_normal( vi ) );
        projectToTangent(fv(vi), vertex_normal,true);
    }
}

void HGF::set_parameters_FES(){
    simulationType=SurfaceFlow;
    
    draw_vertex_velocity=true;
    draw_triangle_velocity=true;
    
    normalize_color = false;
    
    saveMeshColor=false;
    
    th_magnitude=1e-6;
    max_th_nm=Options::doubleValue("max_th_nm") *th_magnitude;
    min_th_nm=Options::doubleValue("min_th_nm")*th_magnitude;
    
    default_th=(max_th_nm + min_th_nm)/2.0;
    
    draw_th_power=1.0;
    draw_vector_scale=0.1;
    
    facePaint=THICKNESS;
    unrendered_region=-1;
    
    triangles_transparent=false;
    
    divergenceMethod=deGoesEtAl;
    
    equation =SoapFilmEquations;
    
    evaporationType=Constant;
    
    preserveLocalThUpdate=Options::boolValue("preserve_local_volumes");;
    preserveLocalThRemesh=false;
    thSmoothingWhenRemeshing=false;
    thSmoothRemeshingCoef=0.5;
    thSmoothEachStep=false;
    thSmoothEachStepCoef=0.99;
    
    maxCrossedEdgesInAdvection = 100;
    
    const std::string vertex_advection_str=Options::strValue("vertex_advection");
    if(vertex_advection_str=="standard"){
        vertexAdvectType=VertStandard;
    }else if(vertex_advection_str=="maccormark"){
        vertexAdvectType=VertMacCormack;
    }else if(vertex_advection_str=="bfecc"){
        vertexAdvectType=VertBFECC;
    }
    
    const std::string triangle_advection_str=Options::strValue("triangle_advection");
    if(triangle_advection_str=="standard"){
        triangleAdvectType=TriStandard;
    }else if(triangle_advection_str=="maccormark"){
        triangleAdvectType=TriMacCormack;
    }else if(triangle_advection_str=="bfecc"){
        triangleAdvectType=TriBFECC;
    }
    
    flow_tension=Options::doubleValue("flow_tension");
    ext_force=Options::doubleValue("ext_force");
    
    vorticityConfinementCoef=Options::doubleValue("vorticity_confinement");
    
    vertical_gravity_scale=Options::doubleValue("vertical_gravity_scale");
    tangential_gravity_scale=Options::doubleValue("tangential_gravity_scale");
    
    thickness_aware_deformation=Options::doubleValue("thickness_aware_deformation");
    deform_standard_thickness=Options::doubleValue("deform_standard_thickness")*th_magnitude;
    if(deform_standard_thickness<0){
        const double ratioToMinTh=1.5;
        deform_standard_thickness=min_th_nm*ratioToMinTh;
    }
    
    externalForceReactionType=VertexLaplacian;
    
    with_evaporation=Options::boolValue("with_evaporation");
    evaporation_per_second = Options::doubleValue("evaporation_per_second");
    
    burstType=Thickness;
    const double burst_ratio = Options::doubleValue("burst_ratio");
    burst_threshold =  min_th_nm+burst_ratio*(max_th_nm-min_th_nm);
    absorbLiquidAftuerBurstCoef=Options::doubleValue("absorb_liquid");
    
    bottom_level=Options::doubleValue("bottom-level");
    
}
#endif //FEOS
