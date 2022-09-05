//
//  HGF_Evolve.cpp
//
//
//  Created by Sadashige Ishida 2017.
//
//

#include <chrono>

#include <igl/cotmatrix.h>
#include <igl/invert_diag.h>
#include <igl/volume.h>
#include <igl/massmatrix.h>
#include <igl/per_vertex_normals.h>

#include "HGF.h"
#include "Options.h"

void HGF::set_intermediate_implicit(double dt, Eigen::MatrixXd &U, Eigen::MatrixXi& F, Eigen::MatrixXd &dUdt, Eigen::MatrixXd &NewU) {
    using namespace Eigen;
    size_t nv = U.rows();
    size_t nt = F.rows();
    
    Eigen::SparseMatrix<double>  L;//cotMatrix
    igl::cotmatrix(U, F, L);
    
    // Recompute just mass matrix on each step
    SparseMatrix<double> M;
    igl::massmatrix(U, F, igl::MASSMATRIX_TYPE_BARYCENTRIC, M);
    
    const auto & S = (M - dt * L);
    
    Eigen::SimplicialLLT<Eigen::SparseMatrix<double > > solver(S);
    
    if (accel) {
        NewU = solver.solve(M * U).eval();
        
        MatrixXd D2UDt2 = surface_tension * (NewU - U) / dt;
        
        if (with_gravity) {
            for (int vi = 0; vi < nv; ++vi) {
                D2UDt2.row(vi) += gravity_vec;
            }
        }
        
        dUdt += D2UDt2 * dt;
        NewU = U + dUdt * dt;
    } else {
        NewU = solver.solve(M * U).eval();
    }
}

// This function currently does not support non-manifold junctions.
void HGF::applyIdealGasPressure(const Eigen::MatrixXd &U,const Eigen::MatrixXi& F,Eigen::MatrixXd& D2UDt2) {
    
    using namespace Eigen;
    
    const size_t nv=mesh().nv();
    VectorXd volumes;
    MatrixXd area_matrix;
    
    volumes_and_areas(U, volumes, area_matrix);
    VectorXd nRTs(num_initialBubbles);
    
    for (int iBubble = 0; iBubble < num_initialBubbles; ++iBubble) {
        nRTs[iBubble] = gas_pressure * initialVolumes[iBubble] * surface_tension / std::pow(3.0 / 4.0 * initialVolumes[iBubble] / M_PI, 1.0 / 3.0);
    }
    
    const VectorXd p_vector = nRTs.array() / volumes.array();
        
    MatrixXd Normals;
    igl::per_vertex_normals(U, F, Normals);
    
    for (size_t vi = 0; vi < nv; ++vi) {
        if (mesh().is_vertex_nonmanifold(vi)) {
            continue;
        }
        
        const size_t tri = mesh().m_vertex_to_triangle_map[vi][0];
        LosTopos::Vec2i v_labels = mesh().get_triangle_label(tri);
        int region0 = (v_labels)[0];
        int region1 = (v_labels)[1];
        
        int non_air_region;
        bool both_non_air = false;
        
        if (region0 <= AIR) {
            //Both sides of the vertex are open regions.
            if (region1 <= AIR) {
                continue;
            }
            
            non_air_region = region1;
        } else if (region1 <= AIR) {
            non_air_region = region0;
        } else {
            both_non_air = true;
        }
        
        if (!both_non_air) {
            D2UDt2.row(vi) += p_vector[non_air_region - region_offset] * Normals.row(vi);
        }
        
        if (both_non_air) {
            int large_region_ind = region0 > region1 ? region0 : region1;
            int small_region_ind = region0 < region1 ? region0 : region1;
            
            D2UDt2.row(vi) += (p_vector[large_region_ind - region_offset] - p_vector[small_region_ind - region_offset]) * Normals.row(vi);
        }
    }
    
}

void HGF::set_intermediate_symplectic_euler(double dt, Eigen::MatrixXd &U, Eigen::MatrixXi& F, Eigen::MatrixXd &dUdt, Eigen::MatrixXd &NewU) {
    using namespace Eigen;
    size_t nv = U.rows();
    size_t nt = F.rows();
    
    MatrixXd dAdx;
    
    Eigen::SparseMatrix<double>  L;//cotMatrix
    
    igl::cotmatrix(U, F, L);
    
    SparseMatrix<double> M, Minv;
    
    igl::massmatrix(U, F, igl::MASSMATRIX_TYPE_BARYCENTRIC, M);
    igl::invert_diag(M, Minv);
    // Laplace-Beltrami of position
    dAdx = -Minv * (L * U);
    
    if (accel) {
        MatrixXd D2UDt2 = -surface_tension * dAdx;
        
        if (with_gravity and vertical_gravity_scale>0.0) {
            
#ifdef FEOS
            MatrixXd projectedGravityMat;
            projectedGravityMat.setZero(nv,3);
            for (int vi = 0; vi < nv; ++vi) {
                
                if(mesh().is_vertex_nonmanifold(vi)){
                    continue;
                }
                projectedGravityMat.row(vi)=vertical_gravity_scale*gravity_vec.dot(vertex_normals[vi])*vertex_normals[vi];
                
            }
            std::vector<size_t>t_junctions,q_junctions;
            nonManifoldJunctions(t_junctions, q_junctions);
            for(const size_t& t_junction:t_junctions){
                
                std::vector<size_t>adjacent_vertices;
                mesh().get_adjacent_vertices(t_junction, adjacent_vertices);
                
                size_t num_valid_adjacent_vertices=0;
                for (size_t ad_v:adjacent_vertices) {
                    if ( not mesh().is_vertex_nonmanifold(ad_v)) {
                        projectedGravityMat.row(t_junction)  += projectedGravityMat.row(ad_v);
                        ++num_valid_adjacent_vertices;
                    }
                }
                if(num_valid_adjacent_vertices>0){
                    projectedGravityMat/=num_valid_adjacent_vertices;
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
                        projectedGravityMat.row(q_junction) += projectedGravityMat.row(ad_v);
                    }
                }
            }
            D2UDt2+=projectedGravityMat;
            
            
#else //FEOS
            for (int vi = 0; vi < nv; ++vi) {
                D2UDt2.row(vi) += gravity_vec;
            }
#endif //FEOS
            
        }
        
#ifdef FEOS
        if(thickness_aware_deformation>0.0){
            VectorXd Th=thvv();
            size_t num_rv=real_vertices.size();
            for (size_t rvi = 0; rvi < num_rv; ++rvi) {
                const size_t vi=real_vertices[rvi];
                
                double deformMultiplier=deform_standard_thickness/Th[vi];
                if(thickness_aware_deformation<1.0){
                    deformMultiplier=std::pow(deformMultiplier, thickness_aware_deformation);
                }
                D2UDt2.row(vi)*= deformMultiplier;
            }
        }
#endif //FEOS
        
        dUdt += D2UDt2 * dt;
        NewU = U + dUdt * dt;
        
    } else {
        double small_scale = 1.0;
        auto DU = -small_scale * surface_tension * dAdx * dt;
        NewU = U + DU;
    }
}

void HGF::stepHGF(double dt) {
    using namespace Eigen;
    
    size_t nv = m_st->m_mesh.nv();
    size_t nt = m_st->m_mesh.nt();
    MatrixXd U = MatrixXd::Zero(nv, 3);
    MatrixXi F = MatrixXi::Zero(nt, 3);
    MatrixXd dUdt = MatrixXd::Zero(nv, 3);
    MatrixXd NewU;
    
    //Convert vertices and faces to U,F, duDt.
    for (size_t i = 0; i < nv; i++) {
        U.row(i) << m_st->pm_positions[i][0], m_st->pm_positions[i][1], m_st->pm_positions[i][2];
        
        //Do not use m_st->pm_velocities.
        //It's not defraged by LosTopos.
        //Instead use "vel" function.
        dUdt.row(i) = vel(i);
    }
    
    for (size_t i = 0; i < nt; i++) {
        F.row(i) << m_st->m_mesh.m_tris[i][0], m_st->m_mesh.m_tris[i][1], m_st->m_mesh.m_tris[i][2];
    }
    
    auto update = [this, &U, &F, &dUdt, &NewU, dt]() {
        if (implicit_scheme) {
            set_intermediate_implicit(dt, U, F, dUdt, NewU);
        } else {
            set_intermediate_symplectic_euler(dt, U, F, dUdt, NewU);
        }
    };
    
    static int count = 0;
    ++count;
    
    if (logging_detailed_time and count > 1) {
        // general time stepping
        std::chrono::system_clock::time_point start, end;
        start = std::chrono::system_clock::now();
        
        update();
        
        end = std::chrono::system_clock::now();
        double elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count(); //Convert time to ms.
        
        time_dAdx += elapsed;
    } else {
        update();
    }
    
    //Ensure constrained vertices to be fixed.
    //Maybe not necessary.
    int nv_cons = constrainedVertices().size();
    
    for (int cvi = 0; cvi < nv_cons; ++cvi) {
        int cv = m_constrained_vertices[cvi];
        NewU.row(cv) = m_constrained_positions[cvi];
        vel(cv) << 0, 0, 0;
    }
    
    if (logging_detailed_time and count > 1) {
        // general time stepping
        std::chrono::system_clock::time_point start, end;
        start = std::chrono::system_clock::now();
        
        correct_volume(NewU, F);
        
        end = std::chrono::system_clock::now();
        double elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count(); //Convert time to ms.
        
        time_volume_correction += elapsed;
    } else {
        correct_volume(NewU, F);
    }
    
    dUdt = (NewU - U) / dt;
    
    dUdt *= damp;
    
    //Convert U and dUdt to vertices and faces.
    for (size_t i = 0; i < nv; i++) {
        m_st->pm_newpositions[i] = vc(NewU.row(i));
        //Do not use m_st->pm_velocities.
        //It's not defraged by LosTopos.
        //Instead use "vel" function.
        vel(i) = dUdt.row(i);
    }
    
    for (int cvi = 0; cvi < nv_cons; ++cvi) {
        int cv = m_constrained_vertices[cvi];
        m_st->pm_newpositions[cv] = vc(m_constrained_positions[cvi]);
        m_st->pm_positions[cv] = vc(m_constrained_positions[cvi]);
        
        vel(cv) << 0, 0, 0;
    }
}

//Correct volume of each closed region.
//This is an extension of [Muller 2009]"Fast and Robust Tracking of Fluid Surfaces" for multiple regions.
void HGF::correct_volume(Eigen::MatrixXd &targetU, Eigen::MatrixXi& F) {
    using namespace Eigen;
    
    MatrixXd Delta_d_n;
    computeDelta_d_n(targetU, F, Delta_d_n);
    targetU += Delta_d_n;
}

void HGF::computeDelta_d_n(const Eigen::MatrixXd &targetU, const Eigen::MatrixXi& F, Eigen::MatrixXd &Delta_d_n) {
    using namespace Eigen;
    int nv = (int)targetU.rows();
    
    if (num_initialBubbles == 0) {
        Delta_d_n = MatrixXd::Zero(nv, 3);
        return;
    }
    
    //Make a linear system by measureing volumes and areas;
    VectorXd volumes;
    
    VectorXd d_vector;
    
    if (sparse) {
        SparseMatrix<double> area_matrix_sparse;
        volumes_and_areas_sparse(targetU, volumes, area_matrix_sparse);
        Eigen::SimplicialLLT<Eigen::SparseMatrix<double > > solver(area_matrix_sparse);
        VectorXd volume_difference = initialVolumes - volumes;
        //Solve Linear system to find pressures of bubbles.
        d_vector = solver.solve(volume_difference).eval();
    } else {
        MatrixXd area_matrix;
        volumes_and_areas(targetU, volumes, area_matrix);
        VectorXd volume_difference = initialVolumes - volumes;
        //Solve Linear system to find pressures of bubbles.
        d_vector = area_matrix.fullPivLu().solve(volume_difference);
    }
    
    //Collect t/q-junctions.
    std::vector<bool> is_tq_junc(nv);
    
    for (int vi = 0; vi < nv; ++vi) {
        is_tq_junc[vi] = mesh().is_vertex_nonmanifold(vi);
    }
    
    //Label of the incident regions for each vertex.
    std::vector<LosTopos::Vec2i> Vlabels(nv);
    
    for (int vi = 0; vi < nv; ++vi) {
        if (is_tq_junc[vi]) {
            continue;
        }
        
        assert(mesh().m_vertex_to_triangle_map[vi].size() > 0);
        size_t tri = mesh().m_vertex_to_triangle_map[vi][0];
        Vlabels[vi] = mesh().get_triangle_label(tri);
    }
    
    //compute current normals.
    //Orientations of normals are from small-indexed region to large-indexed region.
    MatrixXd Normals;
    igl::per_vertex_normals(targetU, F, Normals);
    Delta_d_n = MatrixXd::Zero(nv, 3);
    
    //For non-tq-junctions,
    for (int vi = 0; vi < nv; ++vi) {
        if (is_tq_junc[vi]) {
            continue;
        }
        
        const auto &v_labels = Vlabels[vi];
        
        int region0 = (v_labels)[0];
        int region1 = (v_labels)[1];
        
        int non_air_region;
        bool both_non_air = false;
        
        if (region0 <= AIR) {
            //Both sides of the vertex are open regions.
            if (region1 <= AIR) {
                continue;
            }
            
            non_air_region = region1;
        } else if (region1 <= AIR) {
            non_air_region = region0;
        } else {
            both_non_air = true;
        }
        
        if (!both_non_air) {
            Delta_d_n.row(vi) = d_vector[non_air_region - region_offset] * Normals.row(vi);
        }
        
        if (both_non_air) {
            int large_region_ind = region0 > region1 ? region0 : region1;
            int small_region_ind = region0 < region1 ? region0 : region1;
            
            Delta_d_n.row(vi) = (d_vector[large_region_ind - region_offset] - d_vector[small_region_ind - region_offset]) * Normals.row(vi);
        }
    }
    
    //Special care for t/q-junctions.
    //Needs to be true when handling junctions accurately.
    if (delta_dn_tq_juctions) {
        //For triple junctions
        for (size_t v = 0; v < nv; ++v) {
            if (!is_tq_junc[v]) {
                continue;
            }
            
            Delta_d_n.row(v) << 0, 0, 0;
            
            std::vector<size_t>adjacent_vertices;
            mesh().get_adjacent_vertices(v, adjacent_vertices);
            
            size_t num_valid_adjacent_vertices=0;
            for (size_t ad_v:adjacent_vertices) {
                if (!mesh().is_vertex_nonmanifold(ad_v)) {
                    Delta_d_n.row(v) += Delta_d_n.row(ad_v);
                    ++num_valid_adjacent_vertices;
                }
            }
            
            if(num_valid_adjacent_vertices>0){
                Delta_d_n.row(v) /= num_valid_adjacent_vertices;
                
            }
        }
        
        //For quad or higher-junctions
        std::vector<size_t>q_junctions;
        
        for (size_t v = 0; v < nv; ++v) {
            if (!is_tq_junc[v]) {
                continue;
            }
            
            auto incident_edges = mesh().m_vertex_to_edge_map[v];
            int num_incident_edges = incident_edges.size();
            
            int t_junc_counter = 0;
            
            for (size_t ei = 0; ei < num_incident_edges; ++ei) {
                size_t edge = incident_edges[ei];
                
                if (mesh().is_edge_nonmanifold(edge)) {
                    ++t_junc_counter;
                }
            }
            
            if (t_junc_counter > 3 and !m_st->vertex_is_all_solid(v)) {
                q_junctions.push_back(v);
            }
        }
        
        int num_qj = q_junctions.size();
        
        for (size_t qi = 0; qi < num_qj; ++qi) {
            size_t q = q_junctions[qi];
            Delta_d_n.row(q) << 0, 0, 0;
            
            std::vector<size_t>adjacent_vertices;
            mesh().get_adjacent_vertices(q, adjacent_vertices);
            
            const bool consider_only_t_juncs = true;
            
            if (!consider_only_t_juncs) {
                //Take the average of all the neighbor vertices.
                for (size_t ad_v:adjacent_vertices) {
                    Delta_d_n.row(q) += Delta_d_n.row(ad_v);
                }
            } else {
                //Take the average of only the t-junc neighborfoods.
                const auto& incident_edges = mesh().m_vertex_to_edge_map[q];
                int num_incident_edges = incident_edges.size();
                
                for (size_t ei = 0; ei < num_incident_edges; ++ei) {
                    size_t edge = incident_edges[ei];
                    
                    if (mesh().is_edge_nonmanifold(edge)) {
                        size_t ad_v = edge_other_vertex(ei, q);
                        Delta_d_n.row(q) += Delta_d_n.row(ad_v);
                    }
                }
            }
            
            Delta_d_n.row(q) /= adjacent_vertices.size();
        }
    }
}
