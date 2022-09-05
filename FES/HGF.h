//
//  HGF.h
//  
//
//  Created by Fang Da on 10/27/14.
//
//  Eddited by Sadashige Ishida 2017.

#ifndef __MultiTracker__HGF__
#define __MultiTracker__HGF__

#include <iostream>
#include <memory>
#include <unordered_map>
#include <Eigen/SparseCore>

#include <igl/massmatrix.h>

#include "eigenheaders.h"
#include "surftrack.h"

class Sim;
class Scenes;

class HGF : public LosTopos::SurfTrack::SolidVerticesCallback, public LosTopos::T1Transition::VelocityFieldCallback, public LosTopos::SurfTrack::MeshEventCallback
{
    friend class Sim;
    friend class Scenes;

public:
    //////////////////////////
    //Simulation parameters.
    //////////////////////////
    
    bool accel=true;//The flow is parabolic if it is false.
    
    double damp=1.0;
    double surface_tension=1.0;
    double gas_pressure;

    double smooth=0.1;
    bool smooth_from_neighbors=false;
    //The velocity of each vertex is smoothed using its neighborhoods,
    //if it is true.
    
    bool save_mesh;
    std::string saveMeshExtension;
    bool saveMeshChange_y_z;
    bool saveMeshWith_normal;
    bool saveMeshColor;

    bool sparse=false;
    
    bool implicit_scheme=false;
    
    //Used for logging information.
    bool logging_geometry=true;
    bool write_geometry=true;
    
    bool logging_time=false;
    bool logging_detailed_time=false;
    double time_dAdx,time_volume_correction,time_los_topos;
    double ave_time_dAdx,ave_time_volume_correction,ave_time_los_topos;

    //////////////////////////
    //Simulation parameters
    //dynamically used for a time step.
    //////////////////////////
    bool with_gravity;
    double vertical_gravity_scale;
    double tangential_gravity_scale;
    bool with_evaporation;
    double gravity_scale;
    bool add_velocity=false;
    
    double bottom_level;
    
    bool blowing_bubble0=false;
    bool absorbing_bubble0=false;
    
    bool move_right=false;
    bool move_left=false;
    bool do_stepScenes=false;
    
    bool auto_burst;

public:
    int num_initialBubbles;
    int num_currentBubbles;
    int num_region;
    Eigen::VectorXd initialVolumes;
    
    Eigen::MatrixXd Constrained_V;//Constrained vertices for rendering purposes.
    Eigen::MatrixXi Constrained_F;//Constrained triangles for rendering purposes.

    std::vector<std::unique_ptr<Eigen::Vector3d>> velocity_per_bubble;
    Eigen::Vector3d gravity_vec;

    HGF(const std::vector<LosTopos::Vec3d> & vs, const std::vector<LosTopos::Vec3st> & fs, const std::vector<LosTopos::Vec2i> & ls, const std::vector<size_t> & constrained_vertices= std::vector<size_t>(),  const std::vector<Vec3d> & constrained_positions = std::vector<Vec3d>(),const int num_initialBubbles=-1,const bool detect_boundary=false);
    
    ~HGF();
    
protected:
    LosTopos::SurfTrack * m_st;
     
    bool delta_dn_tq_juctions;
    //Special care for t/q-junction
    //when determining the amount of correction
    //for volume preservation.
    
    int region_offset=1;//Offset for open regions
    int AIR;//Maximum index of open regions.
    //Basically, AIR=region_offet-1;
    
    // Velocities of vertices
    LosTopos::NonDestructiveTriMesh::VertexData<Vec3d> * m_v;
    
    
    // constrained vertices
    std::vector<size_t> m_constrained_vertices;
    std::vector<Vec3d> m_constrained_positions;

public:
    //////////////////////////
    //Main functions
    //////////////////////////
    double step(double dt);

    void easy_orientation();
    //Set each triangle label (i,j) to be i>j.
    
    void set_parameters();
    
    void detect_boundary_vertices();
    //Currently not used.
    //Incompatible with a mesh structure with ghost vertices and open regions.

protected:
    //////////////////////////
    //For updating the state
    //////////////////////////
    double actual_dt;
    
    void stepHGF(double dt);
    
    void set_intermediate_symplectic_euler(double dt,Eigen::MatrixXd &U,Eigen::MatrixXi& F,Eigen::MatrixXd &dUdt,Eigen::MatrixXd &NewU);
    
    void set_intermediate_implicit(double dt,Eigen::MatrixXd &U,Eigen::MatrixXi& F,Eigen::MatrixXd &dUdt,Eigen::MatrixXd &NewU);
    
    void smoothVelocity_neighborhood();
    
    //Correct volume of each closed region.
    //This is an extension of [Muller 2009]"Fast and Robust Tracking of Fluid Surfaces" for multiple regions.
    void correct_volume(Eigen::MatrixXd &targetU,Eigen::MatrixXi& F);
    void computeDelta_d_n(const Eigen::MatrixXd &targetU, const Eigen::MatrixXi& F,Eigen::MatrixXd &Delta_d_n);
    void nonManifoldJunctions(std::vector<size_t>& t_junctions,std::vector<size_t>& q_junctions);

    void applyIdealGasPressure(const Eigen::MatrixXd &U,const Eigen::MatrixXi& F,Eigen::MatrixXd& D2UDt2);
    
public:
    ////////////////////
    //Utilities
    ////////////////////
    void meshAsMatrices(Eigen::MatrixXd& V, Eigen::MatrixXi& F, bool withoutGhostGeometry=false);

    void writeMesh_FaceLabel_constrainedVertices(const bool with_imaginary_vertices=false);
    void geometric_information(double time);
    void write_mesh(const std::string& output_directory);
    void write_film_mesh(std::string=".");
    void write_constrained_mesh(std::string=".");
    void face_edge_vertex_style(Eigen::MatrixXi &FE,Eigen::MatrixXi &EV,Eigen::MatrixXd &V);
    //Make data structure where each row i of FE indicates the edges of the i-th triangle, and each row j of EV indicates the vertices of the j-th edge.

    void volumes_and_areas_sparse(const Eigen::MatrixXd &targetU,Eigen::VectorXd &volumes,Eigen::SparseMatrix<double>& area_matrix );
    void volumes_and_areas(const Eigen::MatrixXd &targetU,Eigen::VectorXd &volumes, Eigen::MatrixXd& area_matrix );
    void volumes_and_areas(Eigen::VectorXd &volumes, Eigen::MatrixXd& area_matrix );

    void total_area(double time,bool write=true);
    void tq_junc_acurracy(double time, bool write=true);
    
    void compute_vertex_areas(std::vector<double> &v_areas);
    double total_energy();

    void test_update_mesh_via_LosTopos();
    
public:
    const LosTopos::SurfTrack * surfTrack() const { return m_st; }
    LosTopos::SurfTrack * surfTrack()       { return m_st; }
    const LosTopos::NonDestructiveTriMesh & mesh() const { return m_st->m_mesh; }
    LosTopos::NonDestructiveTriMesh & mesh()       { return m_st->m_mesh; }
    
    Vec3d pos(size_t v) const { return vc(m_st->pm_positions[v]); }
    const Vec3d & vel(size_t vi) const { return (*m_v)[vi]; }
    Vec3d & vel(size_t vi)       { return (*m_v)[vi]; }
    std::vector<Vec3d> vel() const {
        int nv=mesh().nv();
        std::vector<Vec3d> vels(nv);
        for (size_t i = 0; i <  nv; i++)
        { vels[i] = vel(i);}
        return vels;
    }
    VecXd    velv() const {
        int nv=mesh().nv();
        VecXd vels = VecXd::Zero(nv* 3);
        for (size_t i = 0; i < nv; i++)
        { vels.segment<3>(i * 3) = vel(i); }
        return vels;
    }

    const std::vector<size_t> & constrainedVertices() const { return m_constrained_vertices; }
    std::vector<size_t> & constrainedVertices()       { return m_constrained_vertices; }
    const std::vector<Vec3d> & constrainedPositions() const { return m_constrained_positions; }
    std::vector<Vec3d> & constrainedPositions()       { return m_constrained_positions; }

protected:
    //////////////////////////
    //For acquiring information on mesh connectivity.
    //////////////////////////
    
    // convenient mesh topology query
    size_t edge_other_vertex(size_t e, size_t v) const { LosTopos::Vec2st edge = mesh().m_edges[e]; return edge[0] == v ? edge[1] : edge[0]; }
    // convenient mesh geometry query
    double face_area(size_t f)    const { return m_st->get_triangle_area(f); }
    double vert_area(size_t v)    const { double a = 0; for (size_t i = 0; i < mesh().m_vertex_to_triangle_map[v].size(); i++) a += face_area(mesh().m_vertex_to_triangle_map[v][i]); return a / 3; }
    
    Vec3d  face_outward_normal(size_t f)  const { LosTopos::Vec3st t = mesh().m_tris[f]; Vec3d n = (pos(t[1]) - pos(t[0])).cross(pos(t[2]) - pos(t[0])).normalized(); LosTopos::Vec2i l = mesh().get_triangle_label(f); if (l[0] < l[1]) return -n; else return n; }

    //////////////////////////
    //For mesh update via LosTopos
    //////////////////////////
protected:
    // SurfTrack::SolidVerticesCallback method
    bool            generate_collapsed_position(LosTopos::SurfTrack & st, size_t v0, size_t v1, LosTopos::Vec3d & pos);
    bool            generate_split_position(LosTopos::SurfTrack & st, size_t v0, size_t v1, LosTopos::Vec3d & pos);
    LosTopos::Vec3c generate_collapsed_solid_label(LosTopos::SurfTrack & st, size_t v0, size_t v1, const LosTopos::Vec3c & label0, const LosTopos::Vec3c & label1);
    LosTopos::Vec3c generate_split_solid_label(LosTopos::SurfTrack & st, size_t v0, size_t v1, const LosTopos::Vec3c & label0, const LosTopos::Vec3c & label1);
    bool            generate_edge_popped_positions(LosTopos::SurfTrack & st, size_t oldv, const LosTopos::Vec2i & cut, LosTopos::Vec3d & pos_upper, LosTopos::Vec3d & pos_lower);
    bool            generate_vertex_popped_positions(LosTopos::SurfTrack & st, size_t oldv, int A, int B, LosTopos::Vec3d & pos_a, LosTopos::Vec3d & pos_b);
    bool            solid_edge_is_feature(const LosTopos::SurfTrack & st, size_t e);
    
    // T1Transition::VelocityFieldCallback methods
    LosTopos::Vec3d sampleVelocity(LosTopos::Vec3d & pos);
    bool sampleDirectionalDivergence(const LosTopos::Vec3d & pos, const LosTopos::Vec3d & dir, double & output);
    
    // SurfTrack::MeshEventCallback
    void pre_collapse(const LosTopos::SurfTrack & st, size_t e, void ** data);
    void post_collapse(const LosTopos::SurfTrack & st, size_t e, size_t merged_vertex, void * data);
    
    void pre_split(const LosTopos::SurfTrack & st, size_t e, void ** data);
    void post_split(const LosTopos::SurfTrack & st, size_t e, size_t new_vertex, void * data);
    
    void pre_flip(const LosTopos::SurfTrack & st, size_t e, void ** data);
    void post_flip(const LosTopos::SurfTrack & st, size_t e, void * data);
    
    void pre_t1(const LosTopos::SurfTrack & st, size_t v, void ** data);
    void post_t1(const LosTopos::SurfTrack & st, size_t v, size_t a, size_t b, void * data);
    
    void pre_facesplit(const LosTopos::SurfTrack & st, size_t f, void ** data);
    void post_facesplit(const LosTopos::SurfTrack & st, size_t f, size_t new_vertex, void * data);
    
    void pre_snap(const LosTopos::SurfTrack & st, size_t v0, size_t v1, void ** data);
    void post_snap(const LosTopos::SurfTrack & st, size_t v_kept, size_t v_deleted, void * data);
    
    void pre_smoothing(const LosTopos::SurfTrack & st, void ** data);
    void post_smoothing(const LosTopos::SurfTrack & st, void * data);
    
    std::ostream & log() { static std::stringstream ss; return ss; }

    //For FES
    
#ifdef FEOS
public:
    
    double flow_tension;
    double ext_force;
    double evaporation_per_second;
        
    double thickness_aware_deformation;
    double deform_standard_thickness;
    
    double vorticityConfinementCoef;
    
    bool preserveLocalThUpdate;
    bool preserveLocalThRemesh;
    bool thSmoothingWhenRemeshing;
    double thSmoothRemeshingCoef;
    bool thSmoothEachStep;
    double thSmoothEachStepCoef;
    
    bool velocity_dispersion;
    
    double tri_energy;
    double edge_energy;
        
    double th_magnitude;
    double max_th_nm;
    double min_th_nm;
    double default_th;
    
    double draw_th_power;
    double draw_vector_scale;
    
    bool normalize_color;
    
    bool draw_vertex_velocity;
    bool draw_triangle_velocity;
    
    bool draw_auxiliaries;
    
    enum SimulationGeometry{
        vertex,
        edge,
        triangle
    };
    SimulationGeometry simulationGeometry;
    
    static void projectToTangent(Eigen::Vector3d& target, const Eigen::Vector3d& reference, bool preserveNorm=false);
    static void projectToTangent(Eigen::Vector3d& target, const Eigen::Vector3d& targetNormal, const Eigen::Vector3d& referenceNormal);
    void projectPointToTriangle(Eigen::Vector3d& point, const size_t tri_ind);
    void project_flow_v_tangent();
    
    Eigen::VectorXd P_vert;
	Eigen::VectorXd GradPn_edge;
    Eigen::MatrixXd GradP_vert;
    Eigen::MatrixXd GradP_tri;
    
	Eigen::VectorXd DivV_tri;
	Eigen::VectorXd DivV_vert;
    Eigen::VectorXd DivTV_vert;
    Eigen::VectorXd dTdt_tri;
    Eigen::MatrixXd GradT_vert;
    Eigen::MatrixXd GradT_tri;
    Eigen::VectorXd LT_vert;
    
    Eigen::VectorXd TestVector_vert;
    
    void set_parameters_FES();
    
    void init_FES();
    void init_FES_preserveVelocities();
    
    enum Equation{
        EulerEquation,
        //SWE,
		//WaveEquation,
        SoapFilmEquations
    }equation;
    const std::vector<std::string> equation_names =
            { "Euler Equation", "SoapFilm Equations"};
    
    enum SimulationType{
        SurfaceFlow,
        GometricFlow,
        SurfaceFlowAndGeometricFlow,
    }simulationType;
    const std::vector<std::string> simulationType_names =
            { "Surface flow", "Gometric flow", "Surface flow and geometric flow"};

    enum DivergenceMethod{
        deGoesEtAl, // de Goes et al. "Vector Field Processing on Triangle Meshes" Eqn (9)
        TongEtAl // Tong et al. "Discrete Multiscale Vector Field Decomposition" Eqn. (5)
    }divergenceMethod;
    
    int maxCrossedEdgesInAdvection;
        
    enum TriangleAdvectType{
        TriStandard,
        TriBFECC,
        TriMacCormack
    };
    
    TriangleAdvectType triangleAdvectType;

    enum VertexAdvectType{
        VertStandard,
        VertBFECC,
        VertMacCormack
    };
    VertexAdvectType vertexAdvectType;
    
    void smoothVertexTh(const size_t& new_vertex, const double& smoothingCoef);//0 ~ 1.
    void smoothVertexThs(const double& smoothingCoef);//0 ~ 1.
    
    void computeWeightedLaplacian(Eigen::SparseMatrix<double>& Laplacian,  const Eigen::VectorXd& triValues=Eigen::VectorXd(), const bool withoutGhost=false);
    
    double step_FES(double dt);
    double step_EulerEquation(double dt);
    double step_SoapFilmEquations(double dt);
    
    double step_EulerEquation_TriangleTh(double dt);
    
    void computeSurfaceTensionForce(const Eigen::VectorXd& Th, Eigen::VectorXd& DUth  );
    
    void set_flow_boundary();
    
    double advect_vertex(double dt, Eigen::VectorXd& AdvectionTarget, VertexAdvectType advectType=VertStandard);

    double advect_vertex_standard(double dt, Eigen::VectorXd& AdvectionTarget);
    double advect_vertex_BFECC(double dt, Eigen::VectorXd& AdvectionTarget);
    double advect_vertex_MacCormack(double dt, Eigen::VectorXd& AdvectionTarget);
    double advect_vertex_core(double dt, Eigen::VectorXd& AdvectionTarget, std::vector<size_t>* vertices_near_tjunctions=nullptr,std::vector<size_t>* anomalous_vertices=nullptr);

    double advect_vertex_viaTriangle(double dt, Eigen::VectorXd& AdvectionTarget);
    
    double advect_triangle_viaFlux(double dt, Eigen::VectorXd& AdvectionTarget, const Eigen::MatrixXd& TriVelMat=Eigen::MatrixXd());
        
    Eigen::Vector3d findNextProjectedPrevPosition(const size_t tri_ind, const size_t adj_tri_ind, const Eigen::Vector3d& projected_prev_position, const Eigen::Vector3d& segsIntersection,  size_t edge_ind=-1);

    double advect_triangleVelocity(double dt, TriangleAdvectType advectType=TriStandard);
    double advect_triangle_BFECC(double dt, const Eigen::MatrixXd& AdvectionSource, Eigen::MatrixXd& AdvectionTarget );
    double advect_triangle_MacCormack(double dt, const Eigen::MatrixXd& AdvectionSource, Eigen::MatrixXd& AdvectionTarget );
    double advect_triangle_standard(double dt, const Eigen::MatrixXd& AdvectionSource, Eigen::MatrixXd& AdvectionTarget);
    
    void collectValidTrianglesAndIncomingWaters(const double& dt,const size_t tri_ind, const size_t tjunc_edge, std::vector<size_t>& valid_adj_triangles, std::vector<double>& incoming_waters);
        
    double project_tri(double dt, Equation equation);
    double project_tri(double dt, const Eigen::VectorXd& Uth= Eigen::VectorXd(),const Eigen::VectorXd& Th = Eigen::VectorXd() );

    void computeWeightedDivergenceOnVertices(const Eigen::MatrixXd& vectorField, Eigen::VectorXd& divVec, const Eigen::VectorXd& scalarField=Eigen::VectorXd());
    
    void computeRotationOnVertices(const Eigen::MatrixXd& sourceMat, Eigen::MatrixXd& rotMat);

    void computePressureOnVertices(double dt, const Eigen::VectorXd& Uth= Eigen::VectorXd(), const Eigen::VectorXd& Th= Eigen::VectorXd());
        
    // Currently, vorticityConfinement does not support soapfilms with boundaries
    void vorticityConfinement(Eigen::MatrixXd& triMat, const double coef);
    
    void th_laplacian_filter(double filter_coef);
    
    void scalarVert2tri(const Eigen::VectorXd& vertValues, Eigen::VectorXd& triValues);
    void scalarVert2triNaive(const Eigen::VectorXd& vertValues, Eigen::VectorXd& triValues);
    void scalarVert2triBFECC(const Eigen::VectorXd& vertexValues, Eigen::VectorXd& triangleValues, bool limiter=true);

    void vectorVert2tri(const Eigen::MatrixXd& vertValues, Eigen::MatrixXd& triValues, bool local_energy_preserve=false);
    void vert2triGrad(const Eigen::VectorXd& vertValues, Eigen::MatrixXd& triGrad);
    void vectorTri2vert(const Eigen::MatrixXd& triValues, Eigen::MatrixXd& vertValues, bool local_energy_preserve=false);
    
    void scalarTri2vert(const Eigen::VectorXd& triangleValues, Eigen::VectorXd& vertexValues);
    void scalarTri2vertNaive(const Eigen::VectorXd& triangleValues, Eigen::VectorXd& vertexValues);
    void scalarTri2vertBFECC(const Eigen::VectorXd& triangleValues, Eigen::VectorXd& vertexValues,bool limiter=true);

    void thEdge2Triangle();
    
    void thTriangle2Edge();
    void thVertex2Edge();
	void thEdge2Vertex( );
    void velocityVertex2Triangle(bool local_energy_preserve=false);
    
    void velocityVertex2TriangleStandard(bool local_energy_preserve=false);
    void velocityVertex2TriangleBFECC(bool local_energy_preserve=false);
    
    
    void velocityTriangle2Vertex(bool local_energy_preserve=false);
    void velocityTriangle2VertexStandard(bool local_energy_preserve=false);
    void velocityTriangle2VertexBFECC(bool local_energy_preserve=false);

    void velocityTriangle2VertexFLIP(Eigen::MatrixXd& OldTriVelocities, Eigen::MatrixXd&OldVertVelocities);
    
    void interpolateScalarFromVertices(size_t tri_ind , const Eigen::Vector3d&position, const Eigen::VectorXd& sourceVec, double& destValue);
    void InterpolateScalarFromVerticesOfMultipleTriangles(const std::vector<double>& weights, const std::vector<Eigen::Vector3d>& positions, const std::vector<size_t>& tri_inds, const Eigen::VectorXd& sourceVec, double& destValue);
    
    
    void interpolateVectorFromVertices(size_t tri_ind , const Eigen::Vector3d&position, const Eigen::MatrixXd& sourceMat, Eigen::Vector3d& destVelocity, bool local_energy_preserve);
    void interpolateVectorFromVerticesOnWeightedPaths(const std::vector<size_t>& tri_ind , const std::vector<Eigen::Vector3d>& position, const std::vector<double>& backtrace_weights, const Eigen::Matrix3d& begin_local_frame, const std::vector<Eigen::Matrix3d>& end_local_frames, const Eigen::MatrixXd& sourceMat, Eigen::Vector3d& destVelocity, bool local_energy_preserve);
    void interpolateVectorFromVertices(size_t tri_ind , const std::array<double,3>& weights, const Eigen::MatrixXd& sourceMat, Eigen::Vector3d& destVelocity, bool local_energy_preserve);
        
    void interpolateVelocityFromTriangleShiYu(size_t tri_ind , const Eigen::Vector3d&position, const Eigen::MatrixXd& sourceMat, Eigen::Vector3d&velocity);
    
    enum WeightType{
        UV,
        BaryCentric
    };
    void vertexWeights(std::array<double,3>& weights, size_t tri_ind, const Eigen::Vector3d&position, WeightType weightType=UV );
        
    int semiLagrangianFromTriangle( double dt, size_t tri_ind, Eigen::Vector3d& previous_velocity, std::vector<Eigen::Vector3d>& previous_positions, std::vector<double>& weights, std::vector<size_t>& triangle_indices, std::vector<Eigen::Matrix3d>& end_local_frames);
    
    void calcVertexVelocity( size_t v_ind, int t_ind, Eigen::Vector3d& vertex_velocity, const Eigen::MatrixXd& sourceMat, std::unordered_map<size_t, Eigen::Vector3d>& modified_ad_positions);

    void compute_one_ring_vertices(size_t vertex_index, std::vector<size_t>& adjacent_vertices);
    //int findTargetTriangleByTracking_ShiYu_vertex(double dt, size_t v_ind,Eigen::Vector3d& prev_position, double magPrevVelocity=-1.0, int ite=0);
    int semiLagrangianFromVertex(double dt, size_t v_ind,std::vector<Eigen::Vector3d>& previous_positions, std::vector<double>& weights, std::vector<size_t>& triangle_indices, double magPrevVelocity=-1.0, int ite=0);
    
    int rayIntersectsTriangle(const Eigen::Vector3d& position0, const Eigen::Vector3d& position1, const Eigen::Vector3d& position2, const Eigen::Vector3d& rayOrigin, Eigen::Vector3d& hit_point,double& u, double& v, double EPSILON = 1e-6);
    int rayIntersectsTriangle(size_t tri_ind,const Eigen::Vector3d& rayOrigin, const Eigen::Vector3d& rayDirection, Eigen::Vector3d& hit_point,double EPSILON = 1e-6);
    int twoSegmentsIntersect( Eigen::Vector3d& result, const Eigen::Vector3d& A, const Eigen::Vector3d& B, const Eigen::Vector3d& C, const Eigen::Vector3d& D);
    bool inside_triangle2d(size_t tri_ind, const Eigen::Vector3d& P,double EPSILON = 1e-6);
    bool inside_triangle3d(size_t tri_ind, const Eigen::Vector3d& P, double EPSILON = 1e-6);

    int backtraceFromEdge_woFrame( const double& dt, const size_t e_ind, const size_t tri_ind, const Eigen::Vector3d& intersection, const Eigen::Vector3d& projected_prev_position, std::vector<Eigen::Vector3d>& previous_positions, std::vector<double>& weights, std::vector<size_t>& triangle_indices, int numCrossedEdges);
    
    int backtraceFromEdge_withFrame(const double& dt, const size_t e_ind, const size_t adj_tri_ind, const Eigen::Vector3d& segsIntersection, const Eigen::Vector3d& next_projected_prev_position, std::vector<Eigen::Vector3d>& previous_positions, std::vector<double>& weights, std::vector<size_t>& triangle_indices, std::vector<Eigen::Matrix3d>& end_local_frames, int numCrossedEdges);
    
    int backtraceFromVertex(const double dt, const size_t v_ind, const size_t tri_ind, const Eigen::Vector3d& projected_prev_velocity, const std::unordered_map <size_t,Eigen::Vector3d>& modified_ad_positions,Eigen::Vector3d& naive_prev_position,std::vector<Eigen::Vector3d>& previous_positions, std::vector<double>& weights, std::vector<size_t>& triangle_indices);

	void setLocalFrame(const Eigen::Vector3d& localX, size_t tri_ind, Eigen::Matrix3d& end_local_frame);
    
    std::vector<Eigen::Vector3d>edge_centers;
    std::vector<Eigen::Vector3d>edge_t_normals;

    std::vector<Eigen::Vector3d>edge_v_normals;

    std::vector<double> edge_lengths;
    std::vector<double> edge_tnormal_velocities;
	std::vector<double> delta_edge_tnormal_velocities;

	// For interpolation
    std::vector<Eigen::Vector3d> vertex_velocities_interpolation_error;

    std::vector<double> edge_thickness;

    std::vector<Eigen::Vector3d>triangle_centers;
    std::vector<Eigen::Vector3d>triangle_normals;
    std::vector<double>triangle_areas;
    std::vector<Eigen::Vector3d>triangle_velocities;
    std::vector<double>triangle_thickness;
    
    std::vector<Eigen::Vector3d>vertex_normals;
    std::vector<double>vertex_areas;
    
	Eigen::SparseMatrix<double>triangle_distances;
    
    std::vector<std::vector<size_t>>one_ring_neighbors;
    std::unordered_map<size_t,bool>mapVerticesIrregular;
    
	Eigen::VectorXd DivV_vertIn;
	Eigen::VectorXd P_vertIn;

	Eigen::SparseMatrix<double> MeshLaplacianIn;
	Eigen::SparseMatrix<double> MeshLaplacianTransposeIn;
	Eigen::SparseMatrix<double> MLTMLIn; //=MeshLaplacianTransposeIn*MeshLaplacianIn;
	Eigen::SimplicialCholesky<Eigen::SparseMatrix<double >>LaplacianSolverIn;

    Eigen::SparseMatrix<double> MeshLaplacian;
    Eigen::SparseMatrix<double> MeshLaplacianTranspose;
    Eigen::SparseMatrix<double> MLTML; //=MeshLaplacianTranspose*MeshLaplacian;
//    Eigen::SparseLU<Eigen::SparseMatrix<double >>LaplacianSolver;
    Eigen::SimplicialCholesky<Eigen::SparseMatrix<double >>LaplacianSolver;
//    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double >>LaplacianSolver;

    std::vector<size_t> real_triangles;
    std::vector<size_t> ghost_triangles;
    std::vector<int> triangles_to_real;

    std::vector<size_t> ghost_vertices;
    std::vector<size_t> real_vertices;
    std::vector<int> verts_to_real;
    
    std::vector<size_t> non_solid_vertices;
    std::vector<int> verts_to_non_solids;
    
    double initialLiquidVolume;
    Eigen::VectorXd initialLiquidVolumes;
    size_t numConnectedComponents;
    std::unordered_map<int, int>regionsToComponents;

    double computeTotalLiquidVolume(const Eigen::VectorXd& Th=Eigen::VectorXd());
    void computeLiquidVolumes(Eigen::VectorXd& liquidVolumes,const Eigen::VectorXd& Th=Eigen::VectorXd());
    void correctLiquidVolumes(const double dt, Eigen::VectorXd& Th, Eigen::VectorXd& Uth);
    
    void correctLocalThickness(double dt);
    void affineTransformVectorsOnTriangles(Eigen::MatrixXd& VectorsOnTriangle,const Eigen::MatrixXd& OldPosisions,const Eigen::MatrixXd& NewPositions, bool preserveMagnitude=false);
    
	// This is going to be 3*num_t x 3 matrix.
	// GradB vector of tri,vi_t(=0,1,2) is 3*tri+vi_t th row.
	Eigen::MatrixXd GradB;
	//std::unordered_map<std::pair<size_t, size_t>, Eigen::Vector3d> GradB;
        
    LosTopos::NonDestructiveTriMesh::VertexData<Vec3d> * flow_v;

    const Vec3d & fv(size_t vi) const { return (*flow_v)[vi]; }
    Vec3d & fv(size_t vi)       { return (*flow_v)[vi]; }
    std::vector<Vec3d> fv() const {
        size_t nv=mesh().nv();
        std::vector<Vec3d> flow_vels(nv);
        for (size_t i = 0; i <  nv; i++)
        { flow_vels[i] = fv(i);}
        return flow_vels;
    }
    VecXd    fvv() const {
        size_t nv=mesh().nv();
        VecXd flow_vels = VecXd::Zero(nv* 3);
        for (size_t i = 0; i < nv; i++)
        { flow_vels.segment<3>(i * 3) = fv(i); }
        return flow_vels;
    }
    std::vector<Vec3d>    fvvec() const {
        size_t nv=mesh().nv();
        std::vector<Vec3d> flow_vels;
        flow_vels.reserve(nv);
        for (size_t i = 0; i < nv; i++)
        { flow_vels.push_back(fv(i)); }
        return flow_vels;
    }

	MatXd   fvmat() const {
		size_t nv = mesh().nv();
		MatXd flow_vels = MatXd::Zero( nv, 3 );
		for( size_t i = 0; i < nv; i++ )
		{
			flow_vels.row( i ) =fv( i ) ;
		}
		return flow_vels;
	}

	void setfv( const MatXd& fvMat ){
		size_t nv = mesh().nv();
		for( size_t i = 0; i < nv; i++ )
		{
			fv( i ) = fvMat.row( i );
		}
	}
    
    LosTopos::NonDestructiveTriMesh::FaceData<Vec3d> * flow_t;
    const Vec3d & ft(size_t ti) const { return (*flow_t)[ti]; }
    Vec3d & ft(size_t ti)       { return (*flow_t)[ti]; }
    std::vector<Vec3d> ft() const {
        size_t nt=mesh().nt();
        std::vector<Vec3d> flow_vels(nt);
        for (size_t i = 0; i <  nt; i++)
        { flow_vels[i] = ft(i);}
        return flow_vels;
    }
    VecXd    ftv() const {
        size_t nt=mesh().nt();
        VecXd flow_vels = VecXd::Zero(nt* 3);
        for (size_t i = 0; i < nt; i++)
        { flow_vels.segment<3>(i * 3) = ft(i); }
        return flow_vels;
    }
    std::vector<Vec3d>    ftvec() const {
        size_t nt=mesh().nt();
        std::vector<Vec3d> flow_vels;
        flow_vels.reserve(nt);
        for (size_t i = 0; i < nt; i++)
        { flow_vels.push_back(ft(i)); }
        return flow_vels;
    }
    
    MatXd   ftmat() const {
        size_t nt = mesh().nt();
        MatXd flow_vels = MatXd::Zero( nt, 3 );
        for( size_t i = 0; i < nt; i++ )
        {
            flow_vels.row( i ) =ft( i ) ;
        }
        return flow_vels;
    }
    void setft( const MatXd& ftMat ){
        size_t nt = mesh().nt();
        for( size_t i = 0; i < nt; i++ )
        {
            ft( i ) = ftMat.row( i );
        }
    }
    
    LosTopos::NonDestructiveTriMesh::VertexData<double> * th_v;
    const double & thv(size_t vi) const { return (*th_v)[vi]; }
    double & thv(size_t vi)       { return (*th_v)[vi]; }
    std::vector<double> thv() const {
        size_t nv=mesh().nv();
        std::vector<double> ths(nv);
        for (size_t i = 0; i <  nv; i++)
        { ths[i] = thv(i);}
        return ths;
    }
    VecXd    thvv() const {
        size_t nv=mesh().nv();
        VecXd ths = VecXd::Zero(nv);
        for (size_t i = 0; i < nv; i++)
        { ths[i]= thv(i); }
        return ths;
    }
    void setThv(const VecXd& ThVec){
        size_t nv=mesh().nv();
        for (size_t i = 0; i < nv; i++)
        { thv(i)=ThVec[i] ; }
    }
    
    LosTopos::NonDestructiveTriMesh::VertexData<double> * uth_v;
    const double & uthv(size_t vi) const { return (*uth_v)[vi]; }
    double & uthv(size_t vi)       { return (*uth_v)[vi]; }
    std::vector<double> uthv() const {
        int nv=mesh().nv();
        std::vector<double> uths(nv);
        for (size_t i = 0; i <  nv; i++)
        { uths[i] = uthv(i);}
        return uths;
    }
    VecXd    uthvv() const {
        int nv=mesh().nv();
        VecXd uths = VecXd::Zero(nv);
        for (size_t i = 0; i < nv; i++)
        { uths[i]= uthv(i); }
        return uths;
    }
    void setUThv(const VecXd& UthVec){
        size_t nv=mesh().nv();
        for (size_t i = 0; i < nv; i++)
        { uthv(i)=UthVec[i] ; }
    }

    void setTriVelocity(const MatXd& TriVelocity){
        size_t nt=mesh().nt();
        for (size_t i = 0; i < nt; i++)
        { triangle_velocities[i]=(Vec3d)TriVelocity.row(i); }
    }
    MatXd getTriVelocity(){
        size_t nt=mesh().nt();
        MatXd TriVelocity(nt,3);
        for (size_t i = 0; i < nt; i++)
        { TriVelocity.row(i)=triangle_velocities[i]; }
        return TriVelocity;
    }
    
    enum FacePaint
    {
        NONE,
        THICKNESS
    } facePaint;
    
    bool triangles_transparent;
    // index for one region to be transparent
    int unrendered_region;
    void setUnrenderedTrianglesAndEdges();
    std::unordered_map<size_t,bool> unrendered_edges;
    std::unordered_map<size_t,bool> unrendered_triangles;
    
    // Trajectory of back traces for debugging
	std::vector<Vec3d>trace;
    // Path lengths of back traces for debugging
    std::vector<double>path_lengths;
    
    std::unordered_map<size_t, Vec3d>proj_points_around_0th_ver;
    std::vector<size_t> adjacent_vertices_around_0th_ver;
    std::unordered_map<size_t, Vec3d>proj_vels_around_0th_ver;
    std::unordered_map<size_t, Vec3d>proj_center_around_0th_ver;
    Vec3d ave_proj_vels_around_0th_ver;
    //Vec3d proj_prev_velocity2d;
    Vec3d proj_prev_position2d;
    
    bool intersectionOnce{false};
    
    static std::vector<double> subtractVector(const std::vector<double>& vectorA, const std::vector<double>& vectorB){
        assert(vectorA.size()==vectorB.size());
        size_t num_element =vectorA.size();
        std::vector<double> vectorResult;
        vectorResult.reserve(num_element);
        
        for(size_t i=0;i<num_element;++i){
            vectorResult.push_back(vectorA[i]-vectorB[i]);
        }
        return std::move(vectorResult);
    }
    
    static std::vector<Eigen::Vector2d> subtractVector(const std::vector<Eigen::Vector2d>& vectorA, const std::vector<Eigen::Vector2d>& vectorB){
        assert(vectorA.size()==vectorB.size());
        size_t num_element =vectorA.size();
        std::vector<Eigen::Vector2d> vectorResult;
        vectorResult.reserve(num_element);

        for(size_t i=0;i<num_element;++i){
            vectorResult.push_back(vectorA[i]-vectorB[i]);
        }
        return std::move(vectorResult);
    }
    
    static void printVec(const Vec3d& vec,std::string name=""){
        std::cout<<name<<":"<<vec.x()<<","<<vec.y()<<","<<vec.z()<<"\n";
    };
    
    void printEdges(size_t tind){
        auto edges=mesh().m_triangle_to_edge_map[tind];
        std::cout<<edges[0]<<","<<edges[1]<<","<<edges[2]<<"\n";
    };
    
    void printTri(size_t tind,std::string name=""){
        std::cout<<name<<"\n";
        printVec(pos( mesh().m_tris[tind][0]));
        printVec(pos( mesh().m_tris[tind][1]));
        printVec(pos( mesh().m_tris[tind][2]));
        
        printEdges(tind);
        std::cout<<"\n\n";
    };

	static Eigen::Quaterniond RotationBetweenVectors( const Eigen::Vector3d& normalizedStart, const Eigen::Vector3d& normalizedDest );
    static double clamp(const double& value, const double& minValue, const double& maxValue);
    static void clamp(Eigen::VectorXd& Values, const double& minValue, const double& maxValue);
    
	void writeVelocityOfPoint( double time, bool write = true );
    
    enum EvaporationType{
        Constant,
        Proportional
    }evaporationType;
    
    // VertexDivergence and VertexLaplacian yield the same result.
    // TriangleFlux has strong mesh dependency.
    enum ExternalForceReactionType{
        VertexLaplacian,
        VertexDivergence,
        TriangleFlux
    }externalForceReactionType;
    
    void addGravity2Flow( double dt);
    void evaporate(double dt);
    int findConnectedComponents(std::unordered_map<int,int>& outputRegionToConnectedComponentMap);
    
    void colorForRendering(Eigen::MatrixXd& Color,const double power=1.0);
    
    enum BurstType{
        Random,
        Thickness
    }burstType;
    double burst_threshold;
    double absorbLiquidAftuerBurstCoef;
    
    void absorbLiquidFromABubble(size_t targetRegion);
    #endif //FEOS
    
};

#endif /* defined(__MultiTracker__HGF__) */
