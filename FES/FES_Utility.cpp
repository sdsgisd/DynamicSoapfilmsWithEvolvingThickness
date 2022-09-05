//
//  FES_Utility.cpp
//  LosTopos
//
//  Created by Sadashige ISHIDA on 2018/11/18.
//
#include <stdio.h>
#include <array>

#include <Eigen/SparseCore>
#include <igl/invert_diag.h>
#include <igl/massmatrix.h>
#include <igl/grad.h>
#include <igl/cotmatrix.h>
#include <igl/invert_diag.h>
#include <igl/volume.h>
#include <igl/massmatrix.h>
#include <igl/per_vertex_normals.h>
#include <igl/principal_curvature.h>
#include <igl/colormap.h>

#include "HGF.h"

void HGF::computeSurfaceTensionForce(const Eigen::VectorXd& Th, Eigen::VectorXd& DUth  ){
    
    using namespace Eigen;
    
    enum SurfaceTensionForceType{
        LaplacianThickness,
        DoubleLayerMeanCurvatureNormal,
        DoubleLayerPrincipalCurvatures
    }surfaceTensionForceType=LaplacianThickness;
    
    const size_t num_v = mesh().nv();
    
    switch (surfaceTensionForceType) {
            
        case LaplacianThickness:
        {
            const size_t r_numv=real_vertices.size();
            VectorXd Th_in;
            Th_in.resize(r_numv);
            for(size_t rvi=0;rvi<r_numv;++rvi){
                Th_in[rvi]=Th[real_vertices[rvi]];
            }
            VectorXd DUth_in;
            DUth_in =  MeshLaplacianIn * Th_in;
            DUth.setZero(num_v);
            for(size_t rvi = 0;rvi<real_vertices.size();++rvi){
                DUth[real_vertices[rvi]]=DUth_in[rvi];
            }
            break;
        }
        case DoubleLayerMeanCurvatureNormal:
        {
            
            // Make upper layer and lower layer mesh.
            MatrixXd V,VUpper(num_v,3),VLower(num_v,3);
            MatrixXi F;
            meshAsMatrices(V, F);
            MatrixXd N;
            igl::per_vertex_normals(V, F, N);
            // Coef should be 1.0, but due to numerical precision error, we need some exaggeration.
            
            double coef=1.0;1.0/th_magnitude;
            for(size_t vi=0;vi<num_v;++vi){
                VUpper.row(vi)=V.row(vi)+ coef* 0.5*Th[vi]*N.row(vi);
                VLower.row(vi)=V.row(vi)-coef*0.5*Th[vi]*N.row(vi);
            }
            
            // Compute mean curvature normals of upper/lower mesh.
            Eigen::SparseMatrix<double>  LUpper,LLower;//cotMatrix
            SparseMatrix<double> MUpper, MinvUpper, MLower, MinvLower;
            MatrixXd HnUpper,HnLower;
            
            igl::cotmatrix(VUpper, F, LUpper);
            igl::cotmatrix(VLower, F, LLower);
            
            igl::massmatrix(VUpper, F, igl::MASSMATRIX_TYPE_BARYCENTRIC, MUpper);
            igl::massmatrix(VLower, F, igl::MASSMATRIX_TYPE_BARYCENTRIC, MLower);
            
            igl::invert_diag(MUpper, MinvUpper);
            igl::invert_diag(MLower, MinvLower);
            
            // Laplace-Beltrami of position
            HnUpper = MinvUpper * LUpper * VUpper;
            HnLower = MinvLower * LLower * VLower;
            
            DUth.setZero(num_v,3);
            for(size_t vi=0;vi<num_v;++vi){
                DUth[vi]=(HnUpper.row(vi)-HnLower.row(vi)).dot(N.row(vi));
                
            }
            
            break;
        }
        case DoubleLayerPrincipalCurvatures:
        {
            // Make upper layer and lower layer mesh.
            MatrixXd V,VUpper(num_v,3),VLower(num_v,3);
            MatrixXi F;
            meshAsMatrices(V, F);
            MatrixXd N;
            igl::per_vertex_normals(V, F, N);
            double coef=1.0;
            for(size_t vi=0;vi<num_v;++vi){
                VUpper.row(vi)=V.row(vi)+ coef* 0.5*Th[vi]*N.row(vi);
                VLower.row(vi)=V.row(vi)-coef*0.5*Th[vi]*N.row(vi);
            }
            
            // Compute principal curvatures
            MatrixXd PD1Upper,PD2Upper,PD1Lower,PD2Lower;
            VectorXd PV1Upper,PV2Upper,PV1Lower,PV2Lower;
            igl::principal_curvature(VUpper,F,PD1Upper,PD2Upper,PV1Upper,PV2Upper);
            igl::principal_curvature(VLower,F,PD1Lower,PD2Lower,PV1Lower,PV2Lower);
            
            // Take normal component of mean curvature normals.
            // Normal here is the normal of the original surface (V,F).
            VectorXd meanCurvatureUpper=PV1Upper+PV2Upper;
            VectorXd meanCurvatureLower=PV1Lower+PV2Lower;
            MatrixXd NUpper;
            igl::per_vertex_normals(VUpper, F, NUpper);
            MatrixXd NLower;
            igl::per_vertex_normals(VLower, F, NLower);
            
            for(size_t vi=0;vi<num_v;++vi){
                meanCurvatureUpper[vi]*=NUpper.row(vi).dot(N.row(vi));
                meanCurvatureLower[vi]*=NLower.row(vi).dot(N.row(vi));
            }
            DUth = -meanCurvatureUpper+meanCurvatureLower;
            break;
        }
    }
}

void HGF::computeWeightedLaplacian(Eigen::SparseMatrix<double>& Laplacian, const Eigen::VectorXd& triValues, const bool withoutGhost){
    
    using namespace Eigen;
    
    const VectorXd& TriValues=triValues.size()==mesh().nt()? triValues: VectorXd::Constant(mesh().nt(),1.0);
    
    MatrixXd V;
    MatrixXi F;
    
    //Convert vertices and faces to U,F, and T.
    meshAsMatrices(V,F,withoutGhost);
    
    switch(divergenceMethod){
        case deGoesEtAl:{
            
            if(triValues.size()==0){
                // Make cotangent matrix
                igl::cotmatrix( V, F, Laplacian );
                
            }else{
                
                const size_t r_numv=real_vertices.size();
                
                std::vector<Eigen::Triplet<double> > tripletsForLaplacian;
                tripletsForLaplacian.reserve(r_numv*5*3);
                // In above reserve, 5 is approximated num of triangles for each vertex,
                // 3 is num of adjacent vertices for each triangle.
                for(size_t rvi=0;rvi<r_numv;++rvi){
                    const size_t vi=real_vertices[rvi];
                    const auto& tris=mesh().m_vertex_to_triangle_map[vi];
                    for (const auto& tk_v: tris) {
                        
                        if(triangles_to_real[tk_v]==-1){
                            continue;
                        }
                        
                        const auto& tri=mesh().m_tris[tk_v];
                        size_t vi_in_t=-1;
                        for(size_t i=0;i<3;++i){
                            if(tri[i]==vi){
                                vi_in_t=i;
                            }
                        }
                        
                        for(size_t vj_t=0;vj_t<3;++vj_t){
                            
                            const size_t& vj = tri[ vj_t ];
                            const size_t& v_a = tri[ ( vj_t+1 ) % 3 ];
                            const size_t& v_b = tri[ ( vj_t+2 ) % 3 ];
                            
                            size_t e_ind=mesh().get_edge_index(v_a, v_b);
                            const Vector3d edge = pos(v_b)-pos(v_a );
                            Vector3d outerVector=triangle_normals[tk_v].cross(edge);
                            if(outerVector.dot(edge_centers[e_ind]-triangle_centers[tk_v])<0){
                                outerVector*=-1;
                            }
                            tripletsForLaplacian.emplace_back(rvi,verts_to_real[vj], TriValues[tk_v] * GradB.row( tk_v * 3 + vi_in_t).dot(outerVector));
                            
                        }
                    }
                    
                }
                
                Laplacian.resize(r_numv, r_numv); Laplacian.setFromTriplets(tripletsForLaplacian.begin(), tripletsForLaplacian.end());
                Laplacian*=0.5;
                
            }
            
            break;
        }
            
        case TongEtAl:{
            size_t num_v=mesh().nv();
            
            std::vector<Eigen::Triplet<double>> tripletsForLaplacian;
            tripletsForLaplacian.reserve(num_v*5*3);
            // in above reserve, 5 is expected num of triangles for each vertex,
            // 3 is num of vertices for each triangle
            for(size_t vi=0;vi<num_v;++vi){
                const auto& tris=mesh().m_vertex_to_triangle_map[vi];
                for (const auto& tk_v: tris) {
                    
                    const auto& tri=mesh().m_tris[tk_v];
                    size_t vi_in_t=-1;
                    for(size_t i=0;i<3;++i){
                        if(tri[i]==vi){
                            vi_in_t=i;
                        }
                    }
                    
                    const double& tri_area=triangle_areas[tk_v];
                    
                    for(size_t vj_t=0;vj_t<3;++vj_t){
                        const size_t vj=tri[vj_t];
                        tripletsForLaplacian.emplace_back(vi,vj, TriValues[tk_v]*tri_area*GradB.row( tk_v * 3 + vi_in_t).dot(GradB.row( tk_v * 3 + vj_t)));
                        
                    }
                }
                
            }
            
            Laplacian.resize(num_v, num_v); Laplacian.setFromTriplets(tripletsForLaplacian.begin(), tripletsForLaplacian.end());
            
            break;
        }
    }
    
    SparseMatrix<double> M, Minv;
    
    // Use MASSMATRIX_TYPE_VORONOI
    igl::massmatrix( V, F, igl::MASSMATRIX_TYPE_VORONOI, M );
    igl::invert_diag( M, Minv );
    // Laplace-Beltrami of position
    Laplacian=Minv*Laplacian;
    
}

//http://www.opengl-tutorial.org/jp/intermediate-tutorials/tutorial-17-quaternions/#%E4%BA%8C%E3%81%A4%E3%81%AE%E3%83%99%E3%82%AF%E3%83%88%E3%83%AB%E9%96%93%E3%81%AE%E5%9B%9E%E8%BB%A2%E3%81%AE%E8%A6%8B%E3%81%A4%E3%81%91%E6%96%B9
Eigen::Quaterniond HGF::RotationBetweenVectors( const Eigen::Vector3d& normalizedStart, const Eigen::Vector3d& normalizedDest ) {
    using namespace Eigen;
    
    const double eps_for_unit_check = 1e-8;
    assert( std::abs( normalizedStart.norm() - 1.0 ) < eps_for_unit_check );
    assert( std::abs( normalizedDest.norm() - 1.0 ) < eps_for_unit_check );
    
    const double cosTheta = normalizedStart.dot( normalizedDest );
    
    const double eps_for_unopposite_check = 1e-6;
    assert( !( cosTheta < -1.0 + eps_for_unopposite_check ) );
    
    Vector3d rotationAxis = normalizedStart.cross( normalizedDest );
    
    const double s = sqrt( ( 1 + cosTheta ) * 2 );
    const double invs = 1 / s;
    rotationAxis*=invs;
    
    return Quaterniond( s * 0.5, rotationAxis[0], rotationAxis[1], rotationAxis[2] );
    
}

void HGF::projectPointToTriangle(Eigen::Vector3d& point, const size_t tri_ind){
    Eigen::Vector3d hit_point;
    rayIntersectsTriangle( tri_ind, point, triangle_normals[tri_ind], point );
}

void HGF::projectToTangent(Eigen::Vector3d& target, const Eigen::Vector3d& referenceNormal, bool preserveNorm){
    using namespace Eigen;
    
    if (preserveNorm){
        const double norm=target.norm();
        target -= target.dot( referenceNormal ) / referenceNormal.squaredNorm()*referenceNormal;
        if(!target.isZero()){
            target.normalize();
            target*=norm;
        }
    }else{
        target -= target.dot( referenceNormal ) / referenceNormal.squaredNorm()*referenceNormal;
    }
}

void HGF::projectToTangent(Eigen::Vector3d& target, const Eigen::Vector3d& targetNormal, const Eigen::Vector3d& referenceNormal){
    
    // Theoretically they are same, and preserve the length.
    // I numerically verified they are equivalent.
    enum ProjectType{
        Rotation,
        ParallelTransport
    };
    ProjectType projectType=ParallelTransport;
    
    switch (projectType) {
        case Rotation:
            target.noalias()=RotationBetweenVectors(targetNormal, referenceNormal)*target;
            break;
        case ParallelTransport:
            //https://www.cse.wustl.edu/~taoju/zoum/projects/Livewire/Livewire_PG14_Supplemental.pdf
            //target=a, sumNormal = b+b'
            const Eigen::Vector3d sumNormal=targetNormal+referenceNormal;
            
            target.noalias()-= 2.0*(sumNormal.dot(target))* sumNormal/sumNormal.squaredNorm();
            
            break;
    }
    
}

void HGF::vectorTri2vert(const Eigen::MatrixXd& triValues, Eigen::MatrixXd& vertValues, bool local_energy_preserve){
    
    using namespace Eigen;
    const size_t  num_v = mesh().nv();
    vertValues.setZero(num_v, 3);
    Vector3d vertVector;
    std::unordered_map<size_t, Vector3d> modified_ad_positions;
    
    enum Tri2VertType{
        Standard,
        ShiAndYu
    }tri2vertType=Standard;
    
    for (size_t vi = 0; vi < num_v; ++vi) {
        const bool isTjunction=mesh().is_vertex_nonmanifold(vi);
        
        const auto& inc_tris = mesh().m_vertex_to_triangle_map[ vi ];
        
        double energy_around_vertex = 0.0;
        vertVector.setZero();
        
        if(tri2vertType==ShiAndYu and not isTjunction){
            // FIXME: This would crash with boundary.
            
            
            // pick up one triangle
            const size_t ti_v=inc_tris[0];
            calcVertexVelocity(vi, ti_v, vertVector,triValues, modified_ad_positions);
            projectToTangent(vertVector, triangle_normals[ti_v], vertex_normals[vi]);
            vertValues.row(vi)=vertVector;
            continue;
        }
        
        double validVertexArea=0.0;
        for (size_t inc_tri : inc_tris) {
            
            if(triangles_to_real[inc_tri]==-1){
                continue;
            }
            
            //project triangle velocities to vertex plane
            // i.e. subtract normal component
            
            Vector3d projected_tri_vector=triValues.row(inc_tri);
            
            if(not isTjunction){
                
                projectToTangent(projected_tri_vector, triangle_normals[inc_tri], vertex_normals[vi]);
            }
            
            vertVector += triangle_areas[ inc_tri ] / 3.0*projected_tri_vector;
            
            validVertexArea+=triangle_areas[ inc_tri ] / 3.0;
            //vertex_velocity+=projected_tri_velocity;
            
            if (not isTjunction and local_energy_preserve) {
                energy_around_vertex += triValues.row(inc_tri).squaredNorm()*triangle_areas[ inc_tri ] / 3.0;
            }
            
        }
        
        vertVector /= validVertexArea;
        //vertex_velocity/=inc_tris.size();
        
        if (not isTjunction and local_energy_preserve) {
            if(!vertVector.isZero()){
                vertVector.normalize();
                vertVector *= std::pow( energy_around_vertex / validVertexArea, 0.5 );
            }
            
        }
        vertValues.row(vi)=vertVector;
    }
    
}

void HGF::scalarVert2tri(const Eigen::VectorXd &vertValues, Eigen::VectorXd &triValues){
    scalarVert2triNaive(vertValues,triValues);
}

void HGF::scalarVert2triNaive(const Eigen::VectorXd &vertValues, Eigen::VectorXd &triValues){
    
    size_t nt=mesh().nt();
    triValues.setZero(nt);
    for(size_t ti=0;ti<nt;++ti){
        const auto& tri=mesh().m_tris[ti];
        triValues[ti]=1/3.0*(vertValues[tri[0]]+vertValues[tri[1]]+vertValues[tri[2]]);
    }
    
}

void HGF::scalarVert2triBFECC(const Eigen::VectorXd &vertexValues, Eigen::VectorXd &triangleValues, bool limiter){
    
    using namespace Eigen;
    const size_t num_t = mesh().nt();
    
    const double LargeValue=1e8;
    VectorXd MaxNeighborValues=VectorXd::Constant(num_t,-LargeValue);
    VectorXd MinNeighborValues=VectorXd::Constant(num_t,LargeValue);
    if(limiter){
        // Record max and min adjacent verex values for each triangle.
        for(size_t ti=0; ti<num_t;++ti){
            const auto&tri=mesh().m_tris[ti];
            for(size_t vi_t =0;vi_t <3;++vi_t){
                MaxNeighborValues[ti]=std::max(MaxNeighborValues[ti],vertexValues[tri[vi_t]]);
                MinNeighborValues[ti]=std::min(MinNeighborValues[ti],vertexValues[tri[vi_t]]);
            }
        }
    }
    
    const VectorXd&initialVertexValues=vertexValues;
    VectorXd currentVertexValues=initialVertexValues;
    scalarVert2triNaive(initialVertexValues,triangleValues);
    scalarTri2vertNaive(triangleValues,currentVertexValues);
    
    const VectorXd errorVertexValues=currentVertexValues-initialVertexValues;
    scalarVert2triNaive(initialVertexValues-0.5*errorVertexValues,triangleValues);
    
    if(limiter){
        for(size_t ti=0; ti<num_t;++ti){
            clamp(triangleValues[ti],MinNeighborValues[ti],MaxNeighborValues[ti]);
        }
    }
}

void HGF::vectorVert2tri(const Eigen::MatrixXd& vertValues, Eigen::MatrixXd& triValues, bool local_energy_preserve){
    size_t nt=mesh().nt();
    triValues.setZero(nt,3);
    const std::array<double, 3> weights{1.0/3.0,1.0/3.0,1.0/3.0};
    for(size_t ti=0;ti<nt;++ti){
        Eigen::Vector3d destVector;
        interpolateVectorFromVertices(ti, weights, vertValues,destVector, local_energy_preserve);
        triValues.row(ti)=destVector;
    }
}

void HGF::scalarTri2vert(const Eigen::VectorXd& triangleValues, Eigen::VectorXd& vertexValues){
    scalarTri2vertNaive(triangleValues, vertexValues);
}

void HGF::scalarTri2vertNaive(const Eigen::VectorXd& triangleValues, Eigen::VectorXd& vertexValues){
    // take weighted average of incident triangles
    const size_t num_v=mesh().nv();
    vertexValues.setZero(num_v);
    for(size_t vi=0;vi<num_v;++vi){
        const auto& inc_tris=mesh().m_vertex_to_triangle_map[vi];
        double triple_vertex_area=0.0;
        for(const auto& inc_tri:inc_tris){
            vertexValues[vi]+=triangle_areas[inc_tri]*triangleValues[inc_tri];
            triple_vertex_area+=triangle_areas[inc_tri];
        }
        vertexValues[vi]/=triple_vertex_area;
    }
}

void HGF::scalarTri2vertBFECC(const Eigen::VectorXd& triangleValues, Eigen::VectorXd& vertexValues, bool limiter){
    
    using namespace Eigen;
    
    const size_t num_t=mesh().nt();
    const size_t num_v=mesh().nv();
    
    const double LargeValue=1e8;
    VectorXd MaxNeighborValues=VectorXd::Constant(num_v,-LargeValue);
    VectorXd MinNeighborValues=VectorXd::Constant(num_v,LargeValue);
    if(limiter){
        // Record max and min adjacent verex values for each triangle.
        for(size_t ti=0; ti<num_t;++ti){
            const auto&tri=mesh().m_tris[ti];
            for(size_t vi_t =0;vi_t <3;++vi_t){
                MaxNeighborValues[tri[vi_t]]=std::max(MaxNeighborValues[tri[vi_t]],triangleValues[ti]);
                MinNeighborValues[tri[vi_t]]=std::min(MinNeighborValues[tri[vi_t]],triangleValues[ti]);
                
            }
        }
    }
    
    const VectorXd&initialTriangleValues=triangleValues;
    VectorXd currentTriangleValues=initialTriangleValues;
    scalarTri2vertNaive(initialTriangleValues,vertexValues);
    scalarVert2triNaive(vertexValues,currentTriangleValues);
    
    const VectorXd errorTriangleValues=currentTriangleValues-initialTriangleValues;
    scalarTri2vertNaive(initialTriangleValues-0.5*errorTriangleValues,vertexValues);
    
    if(limiter){
        for(size_t vi=0; vi<num_v;++vi){
            clamp(vertexValues[vi],MinNeighborValues[vi],MaxNeighborValues[vi]);
        }
    }
}

void HGF::vert2triGrad(const Eigen::VectorXd& vertValues, Eigen::MatrixXd& triGrad){
    
    using namespace Eigen;
    
    // Basically both yield the same result.
    enum GradType{
        standard,
        igl
    };
    GradType gradType = igl;
    
    switch (gradType) {
        case standard:{
            size_t nt=mesh().nt();
            triGrad.setZero(nt,3);
            
            for(size_t ti=0;ti<nt;++ti){
                const auto& tri=mesh().m_tris[ti];
                triGrad.row(ti)=1/3.0*(GradB.row(ti*3+0)* vertValues[tri[0]]+GradB.row(ti*3+1)*vertValues[tri[1]]+GradB.row(ti*3+2)*vertValues[tri[2]]);
                
            }
            break;
        }
        case igl:{
            MatrixXd TempGrad=triGrad;
            
            Eigen::MatrixXd V;
            Eigen::MatrixXi F;
            Eigen::SparseMatrix<double> GradMat;
            
            meshAsMatrices(V, F);
            // GradMat is 3nt*nv matrix
            igl::grad(V, F, GradMat);
            triGrad = Map<const MatrixXd>((GradMat*vertValues).eval().data(),F.rows(),3);
            break;
        }
    }
    
}

double HGF::computeTotalLiquidVolume(const Eigen::VectorXd& Th){
    
    const Eigen::VectorXd& Thickness=Th.size()>0?Th:thvv();
    double liquidVol=0.0;
    
    for(size_t rti=0;rti<real_triangles.size();++rti){
        const size_t& ti = real_triangles[rti];
        const auto& tri=mesh().m_tris[ti];
        for(size_t vi_t=0;vi_t<3;++vi_t){
            liquidVol+=Thickness[tri[vi_t]]/3.0*triangle_areas[ti];
        }
    }
    
    return liquidVol;
}

void HGF::computeLiquidVolumes(Eigen::VectorXd& liquidVolumes,const Eigen::VectorXd& Th){
    const Eigen::VectorXd& Thickness=Th.size()>0?Th:thvv();
    liquidVolumes.setZero();
    
    for(size_t rti=0;rti<real_triangles.size();++rti){
        const size_t& ti = real_triangles[rti];
        const auto& tri=mesh().m_tris[ti];
        
        const int& large_tri_label=mesh().m_triangle_labels[rti][0];
        
        for(size_t vi_t=0;vi_t<3;++vi_t){
            
            if(regionsToComponents[large_tri_label]==-1){continue;} liquidVolumes[regionsToComponents[large_tri_label]]+=Thickness[tri[vi_t]]/3.0*triangle_areas[ti];
        }
    }
}

void HGF::correctLiquidVolumes(const double dt, Eigen::VectorXd& Th, Eigen::VectorXd& Uth){
    
    using namespace Eigen;
    const size_t num_v=mesh().nv();
    const VectorXd ThBeforeCorrection=Th;
    
    VectorXd ComponentWiseLiquidVolumes=VectorXd::Zero(numConnectedComponents);
    computeLiquidVolumes(ComponentWiseLiquidVolumes,Th);
    
    const VectorXd DeltaComponentWiseLiquidVolumes=initialLiquidVolumes-ComponentWiseLiquidVolumes;
    VectorXd ComponentWiseAreas=VectorXd::Zero(numConnectedComponents);
    for(const auto& ti:real_triangles){
        const int& largerIndex=mesh().m_triangle_labels[ti][0];
        if(regionsToComponents[largerIndex]==-1){continue;}
        ComponentWiseAreas[regionsToComponents[largerIndex]]+=triangle_areas[ti];
    }
    
    const VectorXd DeltaComponentTh=DeltaComponentWiseLiquidVolumes.array()/ComponentWiseAreas.array();
    VectorXd DeltaTh=VectorXd::Zero(num_v);
    
    for(const size_t& vi:real_vertices){
        auto const& inc_tris=mesh().m_vertex_to_triangle_map[vi];
        
        // vertex area - area of ghost triangles.
        double valid_vertex_area=0.0;
        for(const auto& ti_v:inc_tris){
            if(triangles_to_real[ti_v]==-1){continue;}
            valid_vertex_area+= triangle_areas[ti_v];
            const int& largerIndex=mesh().m_triangle_labels[ti_v][0];
            DeltaTh[vi]+=  DeltaComponentTh[regionsToComponents[largerIndex]]* triangle_areas[ti_v];
        }
        DeltaTh[vi]/=valid_vertex_area;
        
    }
    Th+=DeltaTh;
    clamp(Th,min_th_nm,max_th_nm);
    // adjust uth according to adjusted th.
    for(size_t vi=0;vi<num_v;++vi){
        Uth[vi]+=(Th[vi]-ThBeforeCorrection[vi])/dt;
    }
    
}


double HGF::clamp(const double& value, const double& minValue, const double& maxValue){
    return std::min(std::max(value,minValue),maxValue);
}

void HGF::affineTransformVectorsOnTriangles(Eigen::MatrixXd& VectorsOnTriangles,const Eigen::MatrixXd& OldPosisions, const Eigen::MatrixXd& NewPositions, bool preserveMagnitude){
    
    using namespace Eigen;
    const size_t num_t=mesh().nt();
    Vector3d dummyVector;
    const double eps1=1e-12,eps2=1e-6;
    
    for(size_t ti=0;ti<num_t;++ti){
        
        if(triangles_to_real[ti]==-1){
            VectorsOnTriangles.row(ti)<<0.0,0.0,0.0;
            continue;
        }
        
        if(VectorsOnTriangles.row(ti).squaredNorm()<eps1){
            continue;
        }
        
        const auto& tri=mesh().m_tris[ti];
        
        const double originalMagnituide=VectorsOnTriangles.row(ti).norm();
        
        // Express triangle vector by coefficients u and v.
        double u=0,v=0;
        const Vector3d& old_pos0 = OldPosisions.row( tri[0] );
        const Vector3d& old_pos1 = OldPosisions.row( tri[1] );
        const Vector3d& old_pos2 = OldPosisions.row( tri[2] );
        rayIntersectsTriangle( old_pos0, old_pos1, old_pos2, old_pos0+Vector3d(VectorsOnTriangles.row(ti)), dummyVector, u, v, eps2);
        
        Vector3d TransformedVector= u*(NewPositions.row(tri[1])-NewPositions.row(tri[0]))+v*(NewPositions.row(tri[2])-NewPositions.row(tri[0]));
        if(preserveMagnitude and !TransformedVector.isZero()){
            TransformedVector.normalize();
            TransformedVector*=originalMagnituide;
        }
        VectorsOnTriangles.row(ti)=TransformedVector;
        
    }
    
}

void HGF::correctLocalThickness(double dt){
    using namespace Eigen;
    // Compute new vertex areas;
    size_t num_v=mesh().nv();
    size_t num_t=mesh().nt();
    
    const VectorXd OldTH=thvv();
    
    enum CorrectionMethod{
        AreaBased,
        FluxBased,
    }correctionMethod=AreaBased;
    
    Eigen::MatrixXi F;
    Eigen::MatrixXd OldV;
    meshAsMatrices(OldV, F);
    SparseMatrix<double>OldM;
    igl::massmatrix(OldV, F, igl::MassMatrixType::MASSMATRIX_TYPE_VORONOI ,OldM);
    
    switch(correctionMethod){
        case AreaBased:{
            
            Eigen::VectorXd DeltaTh(num_v);
            
            Eigen::MatrixXd NewV(num_v,3);
            for (size_t vi = 0; vi < num_v; vi++)
            {
                NewV.row(vi)=vc(surfTrack()->pm_newpositions[vi]);
            }
            SparseMatrix<double>NewM;
            igl::massmatrix(NewV, F, igl::MassMatrixType::MASSMATRIX_TYPE_VORONOI ,NewM);
            for (size_t vi = 0; vi < num_v; vi++)
            {
                DeltaTh[vi]=OldM.coeffRef(vi, vi)/NewM.coeffRef(vi, vi)*thv(vi) - thv(vi);
                
            }
            
            const int numSmoothing=3;
            const double smoothCoef=0.95;
            for(int iSmooth=0;iSmooth<numSmoothing;++iSmooth){
                const VectorXd DeltaThRaw=DeltaTh;
                
                for( const size_t& vi:real_vertices )
                {
                    std::vector<size_t>adjacent_vertices;
                    mesh().get_adjacent_vertices(vi, adjacent_vertices);
                    size_t numValidAdjacentVertices=0;
                    double neighborAverage=0.0;
                    for(const size_t& adj_v:adjacent_vertices){
                        if(verts_to_real[adj_v]!=-1){
                            ++numValidAdjacentVertices;
                            neighborAverage+=DeltaThRaw[adj_v];
                        }
                    }
                    if(numValidAdjacentVertices>0){
                        neighborAverage/=numValidAdjacentVertices;
                    }
                    DeltaTh[vi]=(1.0-smoothCoef)*DeltaThRaw[vi]+smoothCoef*neighborAverage;
                }
            }
            
            double correctionCoef=1.0;
            setThv(OldTH+correctionCoef*DeltaTh);
            
            break;
        }
        case FluxBased:{
            // Prepare triangle's Th V where V is vertex velocity not flow velocity.
            MatrixXd Vel_vert=MatrixXd(num_v,3);
            for(size_t vi=0;vi<num_v;++vi){
                //Vel_vert.row(vi)=vc(m_st->pm_newpositions[vi]-m_st->pm_positions[vi]);
                Vel_vert.row(vi)=vel(vi);
            }
            
            MatrixXd TVel_tri=MatrixXd::Zero(num_t, 3);
            for(size_t ti=0;ti<num_t;++ti){
                const auto& tri=mesh().m_tris[ti];
                for(size_t vi_t=0;vi_t<3;++vi_t){
                    const size_t& vert=tri[vi_t];
                    TVel_tri.row(ti)+=1/3.0*Vel_vert.row(vert);
                }
                Vector3d projectedTVel=TVel_tri.row(ti);
                
                projectToTangent(projectedTVel, vc(surfTrack()->get_triangle_normal(ti)));
                TVel_tri.row(ti)=projectedTVel;
            }
            VectorXd TriangleT;
            scalarVert2tri( thvv(), TriangleT );
            
            VectorXd DivTVel_vert;
            computeWeightedDivergenceOnVertices( TVel_tri, DivTVel_vert, TriangleT );
            
            // finally times dt
            for(size_t vi=0;vi<num_v;++vi){
                thv(vi)+=DivTVel_vert[vi]*dt;//*OldM.coeffRef(vi, vi);
            }
            
            break;
        }
    }
    
    // Clamph Th
    for(size_t vi=0;vi<num_v;++vi){
        thv(vi)=clamp(thv(vi), min_th_nm, max_th_nm);
    }
    
    // Update uth.
    const bool updateUTh=0;
    if(updateUTh){
        setUThv(uthvv()+(thvv()-OldTH)/dt);
    }
    
}

void HGF::clamp(Eigen::VectorXd& Values, const double& minValue, const double& maxValue){
    size_t numValues=Values.size();
    for(size_t i=0;i<numValues;++i){
        Values[i]=clamp(Values[i],minValue,maxValue);
    }
}

void HGF::smoothVertexTh(const size_t& vertex, const double& smoothingCoef){
    double neighborhood_mean_th=0.0;
    
    std::vector<size_t> adj_verts;
    mesh().get_adjacent_vertices( vertex, adj_verts);
    for (const size_t& vother: adj_verts)
    {
        neighborhood_mean_th += thv(vother);
    }
    if (not adj_verts.empty()){
        neighborhood_mean_th/=adj_verts.size();
    }
    
    thv(vertex) =(1.0-smoothingCoef)*thv(vertex)+smoothingCoef*neighborhood_mean_th;
}

void HGF::smoothVertexThs(const double& smoothingCoef){
    
    for(size_t vi=0, num_v=mesh().nv();vi<num_v;++vi){
        smoothVertexTh(vi, smoothingCoef);
    }
    
}

// Compute the velocity on a vertex (this is normally stored on triangles)
// Method: see [Shi and Yu, Inviscid and incompressible fluid simulation on triangle meshes] Figure5
// Compute v_ind velocity and store at vertex_velocity
// (i) Unfold incident triangles onto 2D plane and compute vertex velocity as weighted sum of triangles' velocities.
// (ii) Retrieve triangle (t_ind)'s 3D velocity from the 2D velocity obtained in (i).
void HGF::calcVertexVelocity( size_t v_ind, int t_ind, Eigen::Vector3d& vertex_velocity, const Eigen::MatrixXd& sourceMat, std::unordered_map <size_t,Eigen::Vector3d>& modified_ad_positions) {
    
    using namespace Eigen;
    const Vector3d vertex_normal = vc( surfTrack()->get_vertex_normal( v_ind ) );
    
    // Get adjacent vertices
    const std::vector<size_t>& adjacent_vertices=one_ring_neighbors[v_ind];
    
    size_t num_ad_v = adjacent_vertices.size();
    
    // Get the sum of angles of adjacentVertices[i], v_ind, adjacentVertices[i+1] for all i.
    double sum_angle = 0;
    for(int i=0; i<adjacent_vertices.size(); i++) {
        Vector3d edge1 = pos(adjacent_vertices[i]) - pos(v_ind);
        Vector3d edge2 = pos(adjacent_vertices[(i+1)%adjacent_vertices.size()]) - pos(v_ind);
        double angle = acos(edge1.normalized().dot(edge2.normalized()));
        sum_angle += angle;
    }
    
    // Unfold 2D positions of adjacent vertices while keeping the proportion of the angles.
    modified_ad_positions.clear();
    
    // We set adjacent_vertices[0] as 0th vertex around the target vertex
    const Vector3d edge1 = pos(adjacent_vertices[0]) - pos(v_ind);
    modified_ad_positions[adjacent_vertices[0]] = Vector3d(edge1.norm(), 0, 0);
    double crnt_angle = 0;
    
    // Compute other vertex positions.
    // May have some numerical error by acos.
    for(int i=1; i<adjacent_vertices.size(); i++) {
        Vector3d edge1 = pos(adjacent_vertices[i-1]) - pos(v_ind);
        Vector3d edge2 = pos(adjacent_vertices[i]) - pos(v_ind);
        double angle = acos(HGF::clamp(edge1.normalized().dot(edge2.normalized()),-1.0,1.0));
        crnt_angle += (2*M_PI/sum_angle)*angle;
        modified_ad_positions[adjacent_vertices[i]] = Vector3d(edge2.norm()*cos(crnt_angle), 0, edge2.norm()*sin(crnt_angle));
    }
    
    // Compute uv of 3D velocity and get 2D velocity using uv.
    std::vector<size_t> adj_verts;
    mesh().get_adjacent_vertices(v_ind, adj_verts);
    
    const auto& inc_tris = mesh().m_vertex_to_triangle_map[ v_ind ];
    int num_inc_tris = inc_tris.size();
    std::unordered_map<size_t, Vector3d> tri_velocities_from_vert2d;
    std::unordered_map<size_t, Vector3d> tri_center2d;
    std::unordered_map<size_t, double> tri_areas2d;
    double vertex_area = 0;
    
    for (int ti = 0; ti < num_inc_tris; ++ti) {
        
        size_t tri_ind = inc_tris[ ti ];
        Vector3d tri_normal = triangle_normals[ tri_ind ];
        const auto& tri = mesh().m_tris[ tri_ind ];
        
        size_t adj_ver_ind[ 2 ];
        for (int vi = 0, i = 0; vi < 3; ++vi) {
            if (tri[ vi ] != v_ind) {
                adj_ver_ind[ i++ ] = tri[ vi ];
            }
        }
        
        const Vector3d& vert_pos0 = pos( v_ind );
        const Vector3d& vert_pos1 = pos( adj_ver_ind[ 0 ] );
        const Vector3d& vert_pos2 = pos( adj_ver_ind[ 1 ] );
        const Vector3d tri_velocity_from_vert3d = sourceMat.size()>0 ? vert_pos0 +  Vector3d(sourceMat.row(tri_ind)): vert_pos0 + triangle_velocities[tri_ind] ;
        
        Vector3d hit_point;
        double u=0, v=0;
        
        rayIntersectsTriangle( vert_pos0, vert_pos1, vert_pos2, tri_velocity_from_vert3d, hit_point, u, v, 1e-6 );
        
        const Vector3d projected_pos1 = modified_ad_positions.at( adj_ver_ind[ 0 ] );
        const Vector3d projected_pos2 = modified_ad_positions.at( adj_ver_ind[ 1 ] );
        
        const Vector3d tri_velocity_from_vert2d = u*projected_pos1 + v*projected_pos2;
        tri_velocities_from_vert2d[tri_ind] = tri_velocity_from_vert2d;
        tri_center2d[tri_ind] = (projected_pos1 + projected_pos2) / 3.0;
        tri_areas2d[tri_ind] = 0.5 * ( projected_pos1.cross( projected_pos2 ) ).norm();
        vertex_area += 1/3.0 * tri_areas2d[tri_ind];
    }
    
    // Compute 2D vertex velocity
    bool local_energy_preserve=true;
    Vector3d ave_tri_velocities2d(0, 0, 0);
    double energy_around_vertex=0;
    for(auto& adj_velAndProj:tri_velocities_from_vert2d){
        size_t tri_ind = adj_velAndProj.first;
        Vector3d tri_velocity_from_vert2d = tri_velocities_from_vert2d[tri_ind];
        ave_tri_velocities2d += tri_areas2d[tri_ind] / 3.0* tri_velocity_from_vert2d;
        
        if (local_energy_preserve) {
            energy_around_vertex += tri_velocity_from_vert2d.squaredNorm()*tri_areas2d[tri_ind] / 3.0;
        }
    }
    
    ave_tri_velocities2d /= vertex_area;
    
    if (local_energy_preserve) {
        if(!ave_tri_velocities2d.isZero()){
            ave_tri_velocities2d.normalize();
            ave_tri_velocities2d *= std::pow( energy_around_vertex / vertex_area, 0.5 );
        }
    }
    
    // Return 2d velocity that is not on a triangle.
    if(t_ind==-1){
        vertex_velocity=ave_tri_velocities2d;
        return;
    }
    
    // Compute 3D vertex velocity
    const auto& tri = mesh().m_tris[ t_ind ];
    
    size_t adj_ver_ind[ 2 ];
    for (int vi = 0, i = 0; vi < 3; ++vi) {
        if (tri[ vi ] != v_ind) {
            adj_ver_ind[ i++ ] = tri[ vi ];
        }
    }
    
    const Vector3d projected_pos1 = modified_ad_positions.at( adj_ver_ind[ 0 ] );
    const Vector3d projected_pos2 = modified_ad_positions.at( adj_ver_ind[ 1 ] );
    Vector3d hit_point;
    double u=0, v=0, epsForRayTrace=1e-6;
    
    rayIntersectsTriangle( Vector3d::Zero(), projected_pos1, projected_pos2, ave_tri_velocities2d, hit_point, u, v, epsForRayTrace);
    
    const Vector3d& vert_pos0 = pos( v_ind );
    const Vector3d& vert_pos1 = pos( adj_ver_ind[ 0 ] );
    const Vector3d& vert_pos2 = pos( adj_ver_ind[ 1 ] );
    
    const Vector3d& e1 = vert_pos1 - vert_pos0;
    const Vector3d& e2 = vert_pos2 - vert_pos0;
    
    vertex_velocity = u*e1 + v*e2;
}

void HGF::vorticityConfinement(Eigen::MatrixXd &triMat, const double coef){
    using namespace Eigen;
    size_t num_v=mesh().nv();
    size_t num_t=mesh().nt();
    
    std::cout<<"Vorticity Confinement\n";
    std::cout<<"Vorticity Confinement Coef:"<<coef<<"\n";
    
    // Compute vorticity matrix
    MatrixXd Vorticity_vert;
    computeRotationOnVertices(triMat, Vorticity_vert);
    
    VectorXd VorticityNorm_vert;
    VorticityNorm_vert.setZero(num_v);
    for(size_t vi:real_vertices){
        VorticityNorm_vert[vi]=Vorticity_vert.row(vi).norm();
    }
    
    MatrixXd GradVorticityNorm_tri;
    vert2triGrad(VorticityNorm_vert,GradVorticityNorm_tri);
    
    MatrixXd N_tri;
    N_tri.setZero(num_t,3);
    for(size_t ti:real_triangles){
        if(not GradVorticityNorm_tri.row(ti).isZero())
        {
            N_tri.row(ti)=GradVorticityNorm_tri.row(ti).normalized();
        }
    }
    
    MatrixXd Vorticity_tri;
    // Probably this must be false to avoid explosion.
    bool local_length_preserve=false;
    vectorVert2tri(Vorticity_vert, Vorticity_tri,local_length_preserve);
    
    for(size_t ti:real_triangles){
        if(std::isnan( Vorticity_tri.coeffRef(ti,0))
           or std::isnan( Vorticity_tri.coeffRef(ti,1))
           or std::isnan( Vorticity_tri.coeffRef(ti,2))
           or std::isnan( N_tri.coeffRef(ti,0))
           or std::isnan( N_tri.coeffRef(ti,1))
           or std::isnan( N_tri.coeffRef(ti,2)) ){
            std::cout<<"Error detected in vortiricy confinement. Ends the function.\n";
            return;
        }
    }
    
    for(size_t ti:real_triangles){
        const Vector3d& N=N_tri.row(ti);
        const Vector3d& Vorticity=Vorticity_tri.row(ti);
        triMat.row(ti) -= coef* N.cross(Vorticity);
    }
    
}

void HGF::evaporate(double dt){
    using namespace Eigen;
    const size_t num_v = mesh().nv();
    
    VectorXd Th=thvv();
    const VectorXd ThicknessBeforeEvaporation=Th;
    
    switch (evaporationType) {
            // remove a constant amount of thickness
        case Constant:{
            const double evaporationAmount=evaporation_per_second*th_magnitude *dt;
            for(const size_t vi:real_vertices){
                Th[vi] -= evaporationAmount;
            }
            break;
        }
            // remove a fixed proportion of thickness
        case Proportional:{
            const double evaporationRate=100.0*evaporation_per_second*th_magnitude*dt;
            for(const size_t vi:real_vertices){
                Th[vi] -= (Th[vi]-min_th_nm) * evaporationRate;
            }
            break;
        }
    }
    clamp(Th, min_th_nm, max_th_nm);
    setThv(Th);
    
    // You may want to modify uthv according to thickness change.
    const bool adaptUth=false;
    if(adaptUth){
        const VectorXd DeltaTh=thvv()-ThicknessBeforeEvaporation;
        for(size_t vi=0;vi<num_v;++vi){
            uthv(vi)+=DeltaTh[vi]/dt;
        }
    }
    
    computeLiquidVolumes(initialLiquidVolumes);
    
}

void HGF::addGravity2Flow( double dt )
{
    using namespace Eigen;
    size_t num_v = mesh().nv();
    size_t num_t = mesh().nt();
    
    const bool air_resist=true;
    const double air_resist_coef=0.01;
    
    // define delta_fv by projecting coef*dt * gravity_vec to the vertex tangent.
    MatrixXd Delta_fv;
    Delta_fv.setZero(num_v, 3);
    for( size_t vi = 0; vi < num_v; ++vi ){
        
        if(mesh().is_vertex_nonmanifold(vi)){
            continue;
        }
        
        Vector3d delta_fv = dt * gravity_vec;
        projectToTangent(delta_fv, vertex_normals[vi],false);
        Delta_fv.row(vi)=delta_fv;
        
    }
    
    // For non-manifold vertices, interpolate delta_fv from neighbors.
    std::vector<size_t>t_junctions,q_junctions;
    nonManifoldJunctions(t_junctions, q_junctions);
    for(const size_t& t_junction:t_junctions){
        
        std::vector<size_t>adjacent_vertices;
        mesh().get_adjacent_vertices(t_junction, adjacent_vertices);
        
        size_t num_valid_adjacent_vertices=0;
        for (size_t ad_v:adjacent_vertices) {
            if (not mesh().is_vertex_nonmanifold(ad_v)) {
                Delta_fv.row( t_junction ) += Delta_fv.row(ad_v);
                
                ++num_valid_adjacent_vertices;
            }
        }
        if(num_valid_adjacent_vertices>0){
            Delta_fv.row( t_junction )/=num_valid_adjacent_vertices;
        }
    }
    for(const size_t& q_junction:q_junctions){
        std::vector<size_t>adjacent_vertices;
        mesh().get_adjacent_vertices(q_junction, adjacent_vertices);
        
        //Take the average of only the t-junc neighborfoods.
        const auto& incident_edges = mesh().m_vertex_to_edge_map[q_junction];
        
        for (const size_t& inc_edge_ind:incident_edges) {
            
            // Consider only t-junctions
            // i.e. ignore contributions from manifold vertices
            if (mesh().is_edge_nonmanifold(inc_edge_ind)) {
                size_t ad_v = edge_other_vertex(inc_edge_ind, q_junction);
                Delta_fv.row(q_junction) += Delta_fv.row(ad_v);
            }
        }
    }
    for( size_t vi = 0; vi < num_v; ++vi ){
        if(verts_to_real[vi]!=-1){
            fv(vi)+=Delta_fv.row(vi);
        }
    }
    
    // define delta_triangle_velocity by projecting coef*dt * gravity_vec to the triangle tangent.
    
    MatrixXd DeltaVel_tri( num_t, 3 );
    for( size_t ti = 0; ti < num_t; ++ti )
    {
        Vector3d delta_TriVelocity = dt *gravity_vec;
        
        projectToTangent(delta_TriVelocity, triangle_normals[ti],false);
        DeltaVel_tri.row( ti ) = delta_TriVelocity;
        triangle_velocities[ti]+= DeltaVel_tri.row( ti );
        
    }
    
    VectorXd DeltaUth;
    VectorXd DeltaUthTemp;
    
    switch (externalForceReactionType) {
        case VertexLaplacian:{
            
            VectorXd GravityPotential(num_v);
            for (size_t vi=0; vi<num_v; ++vi) {
                GravityPotential[vi]=gravity_vec[2]*surfTrack()->pm_positions[vi][2];
            }
            DeltaUth=-tangential_gravity_scale*dt*MeshLaplacian*GravityPotential;
            break;
        }
        case VertexDivergence:{
            
            VectorXd Div_DeltaVel;
            
            computeWeightedDivergenceOnVertices( DeltaVel_tri,Div_DeltaVel );
            DeltaUth = -tangential_gravity_scale*Div_DeltaVel;
            break;
        }
            
        case TriangleFlux:{
            
            VectorXd UTH_EXTforce_tri=VectorXd::Constant(num_t, 1.0);
            advect_triangle_viaFlux(dt, UTH_EXTforce_tri,DeltaVel_tri);
            
            VectorXd UTH_EXTforce_vert;
            
            scalarTri2vertNaive(UTH_EXTforce_tri, UTH_EXTforce_vert);
            DeltaUth=tangential_gravity_scale*(UTH_EXTforce_vert-VectorXd::Constant(num_v, 1.0));
            
            break;
        }
    }
    
    
    // smooth DeltaUth.
    const bool smoothDeltaUthCoef=1.0;
    if(smoothDeltaUthCoef>0.0){
        const VectorXd DeltaUthRaw=DeltaUth;
        for( const size_t& vi:real_vertices )
        {
            std::vector<size_t>adjacent_vertices;
            mesh().get_adjacent_vertices(vi, adjacent_vertices);
            size_t numValidAdjacentVertices=0;
            double neighborAverage=0.0;
            for(const size_t& adj_v:adjacent_vertices){
                if(verts_to_real[adj_v]!=-1){
                    ++numValidAdjacentVertices;
                    neighborAverage+=DeltaUthRaw[adj_v];
                }
            }
            if(numValidAdjacentVertices>0){
                neighborAverage/=numValidAdjacentVertices;
            }
            DeltaUth[vi]=(1.0-smoothDeltaUthCoef)*DeltaUthRaw[vi]+smoothDeltaUthCoef*neighborAverage;
        }
    }
    
    LT_vert=DeltaUth;
    
    for( size_t vi = 0; vi < num_v; ++vi )
    {
        uthv( vi ) += DeltaUth[vi];
        
    }
    
}

void HGF::collectValidTrianglesAndIncomingWaters(const double& dt,const size_t tri_ind, const size_t tjunc_edge, std::vector<size_t>& valid_adj_triangles, std::vector<double>& incoming_water_rates){
    
    using namespace Eigen;
    const auto& adj_tris=mesh().m_edge_to_triangle_map[ tjunc_edge ];
    
    double sum_incoming_water=0.0;
    for(const size_t adj_tri_ind:adj_tris){
        
        if (adj_tri_ind==tri_ind){continue;}
        const double timeDirection = dt>0?1.0:-1.0;
        const Vector3d velocity= timeDirection*triangle_velocities[adj_tri_ind];// //timeDirection*(triangle_areas[adj_tri_ind]*triangle_velocities[adj_tri_ind]+triangle_areas[tri_ind]*triangle_velocities[tri_ind])/(triangle_areas[tri_ind]+triangle_areas[adj_tri_ind]);
        //timeDirection* 0.5*(triangle_velocities[adj_tri_ind]+triangle_velocities[tri_ind]);
        
        const auto& edgeVertices=mesh().m_edges[tjunc_edge];
        Vector3d edge=pos(edgeVertices[1])-pos(edgeVertices[0]);
        Vector3d edgeNormal=(triangle_normals[adj_tri_ind].cross(edge)).normalized();//(triangle_centers[adj_tri_ind]-triangle_centers[tri_ind]).normalized();
        if(edgeNormal.dot(edge_centers[tjunc_edge]-triangle_centers[adj_tri_ind])<0){
            edgeNormal*=-1;
        }
        
        // We do not have to account for edge length since we only need ratios of incoming water and all the incident triangles share the same edge.
        const double incomingWater=velocity.dot(edgeNormal);
        
        if(incomingWater>0){
            valid_adj_triangles.push_back(adj_tri_ind);
            incoming_water_rates.push_back(incomingWater);
            sum_incoming_water+=incomingWater;
        }
    }
    for(size_t valid_adj_triangle=0;valid_adj_triangle<valid_adj_triangles.size();++valid_adj_triangle){
        incoming_water_rates[valid_adj_triangle]/=sum_incoming_water;
    }
}

Eigen::Vector3d HGF::findNextProjectedPrevPosition(const size_t tri_ind, const size_t adj_tri_ind, const Eigen::Vector3d& projected_prev_position,const Eigen::Vector3d&segsIntersection, size_t edge_ind){
    // rotate (projected_prev_position - intersection) with rotation axis edges[e_ind] and
    // amount dihedral angle(edge*(unprojected_prev_position - intersection)plane , tri_ind_tangent plane)
    
    using namespace Eigen;
    
    if(edge_ind==-1){
        edge_ind=mesh().get_common_edge(tri_ind, adj_tri_ind);
    }
    
    VectorXd sourceTriNormal=triangle_normals[ tri_ind ];
    VectorXd destTriNormal=triangle_normals[ adj_tri_ind ];
    
    // Special treatment for non-manifold edge since the orientation of triangle normals is ambiguous.
    if(mesh().is_edge_nonmanifold(edge_ind)){
        
        const Vector3d edge = pos(mesh().m_edges[edge_ind][1])-pos(mesh().m_edges[edge_ind][0]);
        sourceTriNormal=((edge_centers[edge_ind]-triangle_centers[tri_ind]).cross(edge)).normalized();
        destTriNormal=((triangle_centers[adj_tri_ind]-edge_centers[edge_ind]).cross(edge)).normalized();
        
    }
    
    Quaterniond rotation = RotationBetweenVectors( sourceTriNormal, destTriNormal );
    Vector3d next_projected_prev_position = rotation*( projected_prev_position - segsIntersection ) + segsIntersection;
    
    //Slightly modifying the position by projecting next_projected_prev_position to adj_tri's tangent plane to cancel accumulated numerical error by many time execution of this function.
    bool project_to_triangle=false;
    if(project_to_triangle){
        projectPointToTriangle(next_projected_prev_position, adj_tri_ind);
        
        //    Vector3d hit_point;
        //    rayIntersectsTriangle( adj_tri_ind, next_projected_prev_position, triangle_normals[ adj_tri_ind ], hit_point );
        //    next_projected_prev_position = hit_point;
    }
    
    
    return next_projected_prev_position;
}

// given air_index (number up to which the index represents a ghost component) and a list of triangle labels representing neighboring regions,
// return a map of <region,component> pairs, return value is the number of components (not counting the ghost regions)
int HGF::findConnectedComponents(std::unordered_map<int,int>& outputRegionToConnectedComponentMap){
    outputRegionToConnectedComponentMap.clear();
    const std::vector<LosTopos::Vec2i>& regionPairs = mesh().get_triangle_labels();
    std::unordered_map<int,std::set<int> > neighbor_regions;
    
    // for each region (except outer air), construct a set of neighbor regions
    for (auto i = regionPairs.begin(); i != regionPairs.end(); i++) {
        if (i->v[0] > 0) {
            if (neighbor_regions.find(i->v[0]) == neighbor_regions.end()) {
                std::set<int> neighbor_set;
                neighbor_regions.insert({i->v[0],neighbor_set});
            }
            if (i->v[1] > 0) {
                neighbor_regions.at(i->v[0]).insert(i->v[1]);
            }
        }
        if (i->v[1] > 0) {
            if (neighbor_regions.find(i->v[1]) == neighbor_regions.end()) {
                std::set<int> neighbor_set;
                neighbor_regions.insert({i->v[1],neighbor_set});
            }
            if (i->v[0] > 0) {
                neighbor_regions.at(i->v[1]).insert(i->v[0]);
            }
        }
    }
    
    // perform BFS
    std::queue<int> search_queue;
    int current_vertex;
    int current_component = 0;
    bool component_added;
    for (auto i = neighbor_regions.begin(); i != neighbor_regions.end(); i++) {
        search_queue.push(i->first);
        component_added = false;
        while(!search_queue.empty()) {
            current_vertex = search_queue.front();
            search_queue.pop();
            if (outputRegionToConnectedComponentMap.find(current_vertex) == outputRegionToConnectedComponentMap.end()) {
                outputRegionToConnectedComponentMap.insert({current_vertex,current_component});
                component_added = true;
                for (std::set<int>::iterator j = neighbor_regions.at(current_vertex).begin(); j != neighbor_regions.at(current_vertex).end(); j++) {
                    search_queue.push(*j);
                }
            }
        }
        if (component_added) {
            current_component++;
            component_added = false;
        }
    }
    
    // assign outer air region component -1
    outputRegionToConnectedComponentMap.insert({0,-1});
    
    return current_component;
}

void HGF::colorForRendering(Eigen::MatrixXd& Color,const double power){
    using namespace Eigen;
    const size_t num_v=mesh().nv();
    
    switch (facePaint) {
        case  HGF::FacePaint::NONE:{
            Color.setOnes(num_v, 3);
            break;
        }
        case  HGF::FacePaint::THICKNESS:{
            
            Color.resize(num_v,3);
            
            double mid_th,max_th,min_th;
            
            if(normalize_color){
                const VectorXd Thickness=thvv();
                max_th=Thickness.maxCoeff();
                min_th=Thickness.minCoeff();
                mid_th=0.5*(max_th+min_th);
            }else{
                max_th=max_th_nm;
                min_th=min_th_nm;
                mid_th=default_th;
            }
            
            for(size_t vi=0;vi<num_v;++vi){
                const double&raw_th=thv(vi);
                if(raw_th>mid_th){
                    const double normalizedTh = 1.0 - ( raw_th - mid_th ) / ( max_th - mid_th );
                    const double& poweredTh = power==1.0?normalizedTh: std::pow( normalizedTh, power );
                    Color.row(vi)<<poweredTh, poweredTh, 1.0;
                    
                }else{
                    const double normalizedTh = 1.0 - ( mid_th - raw_th ) / ( mid_th - min_th );
                    const double& poweredTh = power==1.0? normalizedTh:std::pow( normalizedTh, power );
                    Color.row(vi)<<1.0, poweredTh, poweredTh;
                    
                }
            }
            
            break;
        }
    }
    
}

void HGF::absorbLiquidFromABubble(size_t targetRegion){
    
    using namespace Eigen;
    const size_t num_v=mesh().nv();
    VectorXd Th=thvv();
    VectorXd ThTri;
    scalarVert2triNaive(Th, ThTri);
    
    // Correct non-manifold vertices that are incident to the region to pop.
    std::set<size_t> boundaryVertices;
    
    double remainingArea=0.0;
    double wholeAreaOfBurstingRegion=0.0;
    double areaOfBoundary=0.0;
    double lostLiquidVolume=0.0;
    
    for (size_t ti : real_triangles)
    {
        LosTopos::Vec2i l = mesh().get_triangle_label(ti);
        if(l[0] <  l[1]){
            std::swap(l[0],l[1]);
        }
        
        if(l[0] == targetRegion || l[1] == targetRegion){
            wholeAreaOfBurstingRegion+=triangle_areas[ti];
        }
        
        if(l[0] == targetRegion and l[1] == AIR){
            // Compute the liquid volume that is lost by popping.
            lostLiquidVolume+=ThTri[ti]*triangle_areas[ti];
            
            for(size_t vi_t=0;vi_t<3;++vi_t){
                const size_t& vert=mesh().m_tris[ti][vi_t];
                if(mesh().is_vertex_nonmanifold(vert)){
                    boundaryVertices.insert(vert);
                    
                }
            }
        }
        
        else if ( (l[0] == targetRegion and l[1] > AIR) or (l[0] > AIR and l[1] == targetRegion) ){
            remainingArea+=triangle_areas[ti];
        }
        
    }
    
    for(const size_t& boundaryVert:boundaryVertices){
        areaOfBoundary+=vertex_areas[boundaryVert];
    }
    
    const double coef=1.0;
    const double deltaVolume= lostLiquidVolume*remainingArea/wholeAreaOfBurstingRegion;
    //        const double deltaVolume= lostLiquidVolume;
    
    const double deltaVolumePerVertex=coef*deltaVolume/areaOfBoundary;
    
    // TODO: do smoothing
    VectorXd DeltaTh;
    DeltaTh.setZero(num_v);
    
    std::cout<<"lostLiquid:"<<lostLiquidVolume<<", gainLiquid::"<<coef*deltaVolume;
    std::cout<<" deltaVolumePerVertex:"<<deltaVolumePerVertex/th_magnitude<<":"<< deltaVolumePerVertex/th_magnitude*wholeAreaOfBurstingRegion/remainingArea<<"nm\n";
    std::cout<<"boundary vertices:";
    
    for(const size_t& boundaryVert:boundaryVertices){
        DeltaTh[boundaryVert]+=deltaVolumePerVertex;
        std::cout<<boundaryVert<<",";
    }
    std::cout<<"\n";
    
    const int numSmoothing=2;
    const double smoothCoef=0.5;
    for(int iSmooth=0;iSmooth<numSmoothing;++iSmooth){
        const VectorXd DeltaThRaw=DeltaTh;
        
        for( const size_t& vi:real_vertices )
        {
            std::vector<size_t>adjacent_vertices;
            mesh().get_adjacent_vertices(vi, adjacent_vertices);
            size_t numValidAdjacentVertices=0;
            double neighborAverage=0.0;
            for(const size_t& adj_v:adjacent_vertices){
                if(verts_to_real[adj_v]!=-1){
                    ++numValidAdjacentVertices;
                    neighborAverage+=DeltaThRaw[adj_v];
                }
            }
            if(numValidAdjacentVertices>0){
                neighborAverage/=numValidAdjacentVertices;
            }
            DeltaTh[vi]=(1.0-smoothCoef)*DeltaThRaw[vi]+smoothCoef*neighborAverage;
        }
    }
    
    Th+=DeltaTh;
    
    HGF::clamp(Th, min_th_nm, max_th_nm);
    setThv(Th);
    
}
