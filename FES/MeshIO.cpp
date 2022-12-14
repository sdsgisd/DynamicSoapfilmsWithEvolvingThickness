//
//  MeshIO.cpp
//
//  Christopher Batty, Fang Da 2014
//
//
//  Modified by Sadashige Ishida 2017.

#include <fstream>
#include "MeshIO.h"
#include <map>

bool MeshIO::save(HGF & hgf, const std::string & filename, bool binary)
{
    LosTopos::SurfTrack & st = *(hgf.surfTrack());
    
    if (binary)
    {
        std::ofstream os(filename.c_str(), std::ios::binary);
        size_t n;
        
        n = hgf.num_region;
        os.write((char *)&n, sizeof (size_t));
        
        n = st.m_mesh.nv();
        os.write((char *)&n, sizeof (size_t));
        for (size_t i = 0; i < n; i++)
        {
            LosTopos::Vec3d x = st.get_position(i);
            os.write((char *)&(x[0]), sizeof (x[0]));
            os.write((char *)&(x[1]), sizeof (x[1]));
            os.write((char *)&(x[2]), sizeof (x[2]));
            
            //TODO: Insert here code to save velocities.
            
        }
        
        n = st.m_mesh.nt();
        os.write((char *)&n, sizeof (size_t));
        for (size_t i = 0; i < n; i++)
        {
            LosTopos::Vec3st t = st.m_mesh.get_triangle(i);
            if (t[0] == t[1])
                continue;
            LosTopos::Vec2i l = st.m_mesh.get_triangle_label(i);
            
            os.write((char *)&(t[0]), sizeof (t[0]));
            os.write((char *)&(t[1]), sizeof (t[1]));
            os.write((char *)&(t[2]), sizeof (t[2]));
            
            os.write((char *)&(l[0]), sizeof (l[0]));
            os.write((char *)&(l[1]), sizeof (l[1]));
        }
        
        n = hgf.constrainedVertices().size();
        os.write((char *)&n, sizeof (size_t));
        for (size_t i = 0; i < n; i++)
        {
            size_t cv = hgf.constrainedVertices()[i];
            os.write((char *)&cv, sizeof (cv));
            
            Vec3d cx = hgf.constrainedPositions()[i];
            os.write((char *)&(cx[0]), sizeof (cx[0]));
            os.write((char *)&(cx[1]), sizeof (cx[1]));
            os.write((char *)&(cx[2]), sizeof (cx[2]));
        }
        
        return os.good();
    } else
    {
        std::ofstream os(filename.c_str());
        os << hgf.num_region << std::endl;  // nregion
        os << st.m_mesh.nv() << std::endl;
        for (size_t i = 0; i < st.m_mesh.nv(); i++)
        {
            os << st.get_position(i);
            
            //TODO: Insert here code to save velocities.
            os << std::endl;
        }
        
        os << st.m_mesh.nt() << std::endl;
        for (size_t i = 0; i < st.m_mesh.nt(); i++)
        {
            LosTopos::Vec3st t = st.m_mesh.get_triangle(i);
            if (t[0] == t[1])
                continue;
            LosTopos::Vec2i l = st.m_mesh.get_triangle_label(i);
            
            os << t << " " << l << std::endl;
        }
        
        os << hgf.constrainedVertices().size() << std::endl;
        for (size_t i = 0; i < hgf.constrainedVertices().size(); i++)
            os << hgf.constrainedVertices()[i] << " " << hgf.constrainedPositions()[i][0] << " " << hgf.constrainedPositions()[i][1] << " " << hgf.constrainedPositions()[i][2] << std::endl;
        
        os.close();
        
        return os.good();
    }
}

bool MeshIO::load(HGF & hgf, const std::string & filename, bool binary)
{
    LosTopos::SurfTrack & st = *(hgf.surfTrack());
    
    std::ifstream test(filename.c_str());
    if (!test.is_open())
    {
        std::cout << "[MeshIO::load] Error: file " << filename << " not found." << std::endl;
        return false;
    }
    
    for (size_t i = 0; i < st.m_mesh.nt(); i++)
    {
        if (st.m_mesh.get_triangle(i)[0] == st.m_mesh.get_triangle(i)[1])
            continue;
        st.remove_triangle(i);
    }
    
    for (size_t i = 0; i < st.m_mesh.nv(); i++)
        st.remove_vertex(i);
    
    if (binary)
    {
        std::ifstream is(filename.c_str(), std::ios::binary);
        
        size_t n;
        is.read((char *)&n, sizeof (size_t));
        size_t nregion = n;
        
        n = st.m_mesh.nv();
        is.read((char *)&n, sizeof (size_t));
        
        st.m_mesh.set_num_vertices(n);
        std::vector<LosTopos::Vec3d> pos(n);
        for (size_t i = 0; i < n; i++)
        {
            LosTopos::Vec3d x;
            is.read((char *)&(x[0]), sizeof (x[0]));
            is.read((char *)&(x[1]), sizeof (x[1]));
            is.read((char *)&(x[2]), sizeof (x[2]));
            pos[i] = x;
            
            //Insert here code to load velocities.
        }
        
        st.m_masses.resize(n);
        for (size_t i = 0; i < n; i++)
            st.m_masses[i] = LosTopos::Vec3d(1, 1, 1);
        
        st.pm_positions = pos;
        st.pm_newpositions = pos;
        st.set_all_remesh_velocities(std::vector<LosTopos::Vec3d>(n, LosTopos::Vec3d(0)));
        
        n = st.m_mesh.nt();
        is.read((char *)&n, sizeof (size_t));
        
        std::vector<LosTopos::Vec3st> tris;
        std::vector<LosTopos::Vec2i> labels;
        for (size_t i = 0; i < n; i++)
        {
            LosTopos::Vec3st t;
            is.read((char *)&(t[0]), sizeof (t[0]));
            is.read((char *)&(t[1]), sizeof (t[1]));
            is.read((char *)&(t[2]), sizeof (t[2]));
            tris.push_back(t);
            
            LosTopos::Vec2i l;
            is.read((char *)&(l[0]), sizeof (l[0]));
            is.read((char *)&(l[1]), sizeof (l[1]));
            labels.push_back(l);
        }
        
        st.m_mesh.replace_all_triangles(tris, labels);
        
        size_t nv = st.m_mesh.m_vertex_to_triangle_map.size();
        st.pm_positions.resize(nv);
        st.pm_newpositions.resize(nv);
        st.pm_velocities.resize(nv);
        st.m_velocities.resize(nv);
        
        st.set_all_positions(pos);
        st.set_all_newpositions(pos);
        
        bool good = is.good();
        
        n = 0;
        is.read((char *)&n, sizeof (size_t));
        if (!is.eof())  // if eof, this means the rec file is older version (not containing constraints info)
        {
            hgf.constrainedVertices().clear();
            hgf.constrainedPositions().clear();
            
            for (size_t i = 0; i < n; i++)
            {
                size_t cv;
                is.read((char *)&cv, sizeof (cv));
                
                Vec3d cx;
                is.read((char *)&(cx[0]), sizeof (cx[0]));
                is.read((char *)&(cx[1]), sizeof (cx[1]));
                is.read((char *)&(cx[2]), sizeof (cx[2]));
                
                hgf.constrainedVertices().push_back(cv);
                hgf.constrainedPositions().push_back(cx);
            }
            
            for (size_t i = 0; i < hgf.constrainedVertices().size(); i++)
                st.m_masses[hgf.constrainedVertices()[i]] = LosTopos::Vec3d(1, 1, 1) * std::numeric_limits<double>::infinity();
            
            good = good && is.good();
        }
        
        is.close();
        
        return good;
    } else
    {
        std::ifstream is(filename.c_str());
        
        size_t n;
        is >> n;
        size_t nregion = n;
        
        is >> n;
        st.m_mesh.set_num_vertices(n);
        std::vector<LosTopos::Vec3d> pos(n);
        for (size_t i = 0; i < n; i++)
        {
            LosTopos::Vec3d x;
            is >> x[0] >> x[1] >> x[2];
            pos[i] = x;
            
            //TODO: Insert here code to load velocities.
        }
        
        st.m_masses.resize(n);
        for (size_t i = 0; i < n; i++)
            st.m_masses[i] = LosTopos::Vec3d(1, 1, 1);
        
        st.pm_positions = pos;
        st.pm_newpositions = pos;
        st.set_all_remesh_velocities(std::vector<LosTopos::Vec3d>(n, LosTopos::Vec3d(0)));
        
        n = st.m_mesh.nt();
        is >> n;
        
        std::vector<LosTopos::Vec3st> tris;
        std::vector<LosTopos::Vec2i> labels;
        for (size_t i = 0; i < n; i++)
        {
            LosTopos::Vec3st t;
            is >> t[0] >> t[1] >> t[2];
            tris.push_back(t);
            
            LosTopos::Vec2i l;
            is >> l[0] >> l[1];
            labels.push_back(l);
        }
        
        st.m_mesh.replace_all_triangles(tris, labels);
        
        size_t nv = st.m_mesh.m_vertex_to_triangle_map.size();
        st.pm_positions.resize(nv);
        st.pm_newpositions.resize(nv);
        st.pm_velocities.resize(nv);
        st.m_velocities.resize(nv);
        
        st.set_all_positions(pos);
        st.set_all_newpositions(pos);
        
        //TODO: Insert here store velocities to hgf.m_v.
        //        for (size_t i = 0; i < vs.mesh().nv(); i++)
        //            vs.vel(i)=velocities[i];
        
        bool good = is.good();
        
        is >> n;
        if (!is.eof())
        {
            hgf.constrainedVertices().clear();
            hgf.constrainedPositions().clear();
            
            for (size_t i = 0; i < n; i++)
            {
                size_t cv;
                is >> cv;
                
                Vec3d cx;
                is >> cx[0] >> cx[1] >> cx[2];
                
                hgf.constrainedVertices().push_back(cv);
                hgf.constrainedPositions().push_back(cx);
            }
            
            for (size_t i = 0; i < hgf.constrainedVertices().size(); i++)
                st.m_masses[hgf.constrainedVertices()[i]] = LosTopos::Vec3d(1, 1, 1) * std::numeric_limits<double>::infinity();
            
            good = good && is.good();
        }
        
        is.close();
        
        return good;
    }
}

bool MeshIO::loadIntoRaw(std::vector<LosTopos::Vec3d> & vs, std::vector<LosTopos::Vec3st> & fs, std::vector<LosTopos::Vec2i> & ls, const std::string & filename, bool binary)
{
    std::ifstream test(filename.c_str());
    if (!test.is_open())
    {
        std::cout << "[MeshIO::load] Error: file " << filename << " not found." << std::endl;
        return false;
    }
    
    if (binary)
    {
        std::ifstream is(filename.c_str(), std::ios::binary);
        
        size_t n;
        is.read((char *)&n, sizeof (size_t));
        
        vs.resize(n);
        for (size_t i = 0; i < n; i++)
        {
            LosTopos::Vec3d x;
            is.read((char *)&(x[0]), sizeof (x[0]));
            is.read((char *)&(x[1]), sizeof (x[1]));
            is.read((char *)&(x[2]), sizeof (x[2]));
            vs[i] = x;
        }
        
        is.read((char *)&n, sizeof (size_t));
        
        fs.clear();
        ls.clear();
        fs.reserve(n);
        ls.reserve(n);
        for (size_t i = 0; i < n; i++)
        {
            LosTopos::Vec3st t;
            is.read((char *)&(t[0]), sizeof (t[0]));
            is.read((char *)&(t[1]), sizeof (t[1]));
            is.read((char *)&(t[2]), sizeof (t[2]));
            fs.push_back(t);
            
            LosTopos::Vec2i l;
            is.read((char *)&(l[0]), sizeof (l[0]));
            is.read((char *)&(l[1]), sizeof (l[1]));
            ls.push_back(l);
        }
        
        is.close();
        
        return is.good();
    } else
    {
        std::ifstream is(filename.c_str());
        
        size_t n;
        is >> n;
        
        vs.resize(n);
        for (size_t i = 0; i < n; i++)
        {
            LosTopos::Vec3d x;
            is >> x[0] >> x[1] >> x[2];
            vs[i] = x;
        }
        
        is >> n;
        
        fs.clear();
        ls.clear();
        fs.reserve(n);
        ls.reserve(n);
        for (size_t i = 0; i < n; i++)
        {
            LosTopos::Vec3st t;
            is >> t[0] >> t[1] >> t[2];
            fs.push_back(t);
            
            LosTopos::Vec2i l;
            is >> l[0] >> l[1];
            ls.push_back(l);
        }
        
        is.close();
        
        return is.good();
    }
}

namespace
{
    struct VertexPatchPair
    {
        size_t v;
        Vec2i rp;   // must satisfy rp[0] < rp[1]
        
        bool operator () (const VertexPatchPair & vpp1, const VertexPatchPair & vpp2) const { return (vpp1.v < vpp2.v || (vpp1.v == vpp2.v && (vpp1.rp[0] < vpp2.rp[0] || (vpp1.rp[0] == vpp2.rp[0] && vpp1.rp[1] < vpp2.rp[1])))); }
    };
    
}

bool MeshIO::saveOBJ(HGF & hgf, const std::string & filename,const bool with_imaginary_vertices, const bool with_normal,const bool change_y_z)
{
    std::cout << "Saving mesh to " << filename << "... ";
    
    // define all the normals (for manifold vertices, average the incident face normals weighted by area; for nonmanifold vertex, one normal instance is created for each manifold patch.)
    std::map<VertexPatchPair, int, VertexPatchPair> vertex_normal_index;
    std::vector<Vec3d> vertex_normals;
    if(with_normal){
        
        for (size_t i = 0; i < hgf.mesh().nv(); i++)
        {
            std::set<Vec2i, Vec2iComp> incident_regions;
            for (size_t j = 0; j < hgf.mesh().m_vertex_to_triangle_map[i].size(); j++)
            {
                LosTopos::Vec2i l = hgf.mesh().get_triangle_label(hgf.mesh().m_vertex_to_triangle_map[i][j]);
                Vec2i rp = (l[0] < l[1] ? Vec2i(l[0], l[1]) : Vec2i(l[1], l[0]));
                incident_regions.insert(rp);
            }
            
            for (std::set<Vec2i, Vec2iComp>::iterator j = incident_regions.begin(); j != incident_regions.end(); j++)
            {
                VertexPatchPair vpp;
                vpp.v = i;
                vpp.rp = *j;
                vertex_normal_index[vpp] = (int)vertex_normals.size();
                Vec3d n(0, 0, 0);
                for (size_t k = 0; k < hgf.mesh().m_vertex_to_triangle_map[i].size(); k++)
                {
                    LosTopos::Vec2i l = hgf.mesh().get_triangle_label(hgf.mesh().m_vertex_to_triangle_map[i][k]);
                    if ((l[0] == (*j)[0] && l[1] == (*j)[1]) || (l[1] == (*j)[0] && l[0] == (*j)[1]))
                    {
                        LosTopos::Vec3st t = hgf.mesh().get_triangle(hgf.mesh().m_vertex_to_triangle_map[i][k]);
                        Vec3d fn = (hgf.pos(t[1]) - hgf.pos(t[0])).cross(hgf.pos(t[2]) - hgf.pos(t[0])) * (l[0] == (*j)[0] ? 1 : -1);
                        n += fn;
                    }
                }
                n.normalize();
                vertex_normals.push_back(n);
            }
        }
        
    }
    
    std::ofstream os(filename.c_str());
    for (size_t i = 0; i < hgf.mesh().nv(); i++){
        
        double x=hgf.pos(i).x(),y=hgf.pos(i).y(),z=hgf.pos(i).z();
        if(change_y_z){
            std::swap(y,z);
        }
        
        os << "v " << x << " " << y<< " " << z<< std::endl;
        
    }
    
    if(with_normal){
        for (size_t i = 0; i < vertex_normals.size(); i++){
        
            double x=vertex_normals[i].x(),y=vertex_normals[i].y(),z=vertex_normals[i].z();
        if(change_y_z){
            std::swap(y,z);
        }
        
        if(hgf.saveMeshColor){
            // Need to flip normal so that meshes are rendered properly with mitsuba.
            x*=-1.0;
            y*=-1.0;
            z*=-1.0;
        }
            
        os << "vn " << x << " " << y << " " << z << std::endl;
            
        }
    }
    
#ifdef FEOS
    
    if(hgf.saveMeshColor){
        Eigen::MatrixXd Color;
        hgf.colorForRendering(Color,1.0);
        const Eigen::MatrixXi IntColor=(256*Color).cast<int>();

        for (size_t i = 0; i < hgf.mesh().nv(); i++){
            os << "vc " << IntColor(i,0)<<" "<<IntColor(i,1)<<" "<<IntColor(i,2) << std::endl;
        }
        
    }else{
        for (size_t i = 0; i < hgf.mesh().nv(); i++){
            os << "th " << hgf.thv( i )/hgf.th_magnitude<< std::endl;
        }
    }
#endif //FEOS
    
    for (size_t i = 0; i < hgf.mesh().nt(); i++)
    {
        
        bool special_case_for_cubicframe=0 ;
        
        if(!with_imaginary_vertices and hgf.surfTrack()->triangle_is_all_solid(i)){
            
            if(!special_case_for_cubicframe){
                continue;
                
            }else{
                LosTopos::Vec3st t = hgf.mesh().get_triangle(i);
                auto abspos0=hgf.pos(t[0]).cwiseAbs();
                auto abspos1=hgf.pos(t[1]).cwiseAbs();
                auto abspos2=hgf.pos(t[2]).cwiseAbs();
                
                if(abspos0.maxCoeff()> 3 or abspos1.maxCoeff()>3 or abspos2.maxCoeff()>3){
                    continue;
                }
                
            }
            
        }
        
        LosTopos::Vec3st t = hgf.mesh().get_triangle(i);
        LosTopos::Vec2i l = hgf.mesh().get_triangle_label(i);
        VertexPatchPair vpp;
        
        if(change_y_z){
            std::swap(t[1],t[2]);
        }
        
        if(with_normal){
            vpp.rp = (l[0] < l[1] ? Vec2i(l[0], l[1]) : Vec2i(l[1], l[0]));
            vpp.v = t[0];   int vn0 = vertex_normal_index[vpp];
            vpp.v = t[1];   int vn1 = vertex_normal_index[vpp];
            vpp.v = t[2];   int vn2 = vertex_normal_index[vpp];
            
            os << "f " << t[0] + 1 << "//" << vn0 + 1 << " " << t[1] + 1 << "//" << vn1 + 1 << " " << t[2] + 1 << "//" << vn2 + 1 << std::endl;
            
        }else{
            os << "f " << t[0] + 1 << " " << t[1] + 1 << " " << t[2] + 1 <<  std::endl;
        }
    }
    
    os.close();
    
    std::cout << "done." << std::endl;
    
    return true;
}

bool MeshIO::savePLY(HGF & hgf, const std::string & filename,const bool with_imaginary_vertices, const bool change_y_z){
    
    std::cout << "Saving mesh to " << filename << "... ";
    
    const size_t numTriangles=with_imaginary_vertices?hgf.mesh().nt():hgf.real_triangles.size();

    std::ofstream os(filename.c_str());
    
    os << "ply" << std::endl;
    os << "format ascii 1.0" << std::endl;
    os << "comment VCGLIB generated" << std::endl;
    os << "element vertex " << hgf.mesh().nv() <<std::endl;
    os << "property float x" << std::endl;
    os << "property float y" << std::endl;
    os << "property float z" << std::endl;
    os << "property uchar red" << std::endl;
    os << "property uchar green" << std::endl;
    os << "property uchar blue" << std::endl;
    os << "property uchar alpha" << std::endl;
    os << "element face " << numTriangles << std::endl;
    os << "property list uchar int vertex_indices" << std::endl;
    os << "end_header" << std::endl;
    
    for (size_t i = 0; i < hgf.mesh().nv(); i++){
        double x=hgf.pos(i).x(),y=hgf.pos(i).y(),z=hgf.pos(i).z();
        if(change_y_z){
            std::swap(y,z);
        }
        
        double raw_th=hgf.thv(i);
        if(1){
            raw_th = 255*(raw_th - hgf.min_th_nm)/(hgf.max_th_nm - hgf.min_th_nm);
//            std::cout<<raw_th<<","<<int(255*(raw_th - hgf.min_th_nm)/(hgf.max_th_nm - hgf.min_th_nm))<<"\n";
            
            raw_th=HGF::clamp(raw_th, 0,255);
        }
        
        //            os << x << " " << y << " " << z << " " << 0*std::round(raw_th) << " " << 0 << " " << 0 << " " << 255 << std::endl;
        
        os << x << " " << y << " " << z << " " << int(raw_th) << " " << 255 << " " << 255 << " " << 255 << std::endl;
        
    }
    
    for (size_t i = 0; i < hgf.mesh().nt(); i++)
    {
        if(!with_imaginary_vertices and hgf.surfTrack()->triangle_is_all_solid(i)){continue;}
        LosTopos::Vec3st t = hgf.mesh().get_triangle(i);
        if(change_y_z){
            std::swap(t[1],t[2]);
        }
        
        //ply file does not require "+1" unlike obj.
        os << "3 " << t[0] << " " << t[1] << " " << t[2] <<  std::endl;
    }
    
    os.close();
    
    std::cout << "done." << std::endl;
    
    return true;
    
}
