//
//  Sim.cpp
//
//  Fang Da 2014
//
//  Editted by Sadashige Ishida 2017

#include <sstream>
#include <iomanip>
#include <chrono>
#include "Sim.h"
#include "Options.h"
#include "MeshIO.h"
#include <cmath>
#include "eigenheaders.h"

#ifndef WIN32
#include <sys/time.h>
#include <unistd.h>
#else
#include <direct.h>
#include <windows.h>
#endif //WIN32

#include <sys/types.h>
#include <sys/stat.h>
#include "YImage.h"
#include "Colormap.h"
#include "PRRenderer.h"

#ifdef __APPLE__
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#include <GLUT/glut.h>
#elif __linux__
#include <GL/glew.h>
#include <GL/glut.h>
#include <GL/gl.h>
#include <GL/glu.h>
#elif WIN32
#include <windows.h>
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>
#endif

Sim::Sim(bool verbose) :
m_verbose(verbose),
m_scene("unspecified"),
m_output_directory(""),
hgf(NULL),
m_dt(0),
m_time(0),
m_frameid(0),
m_finished(false),
m_nearest_vertex(-1),
m_nearest_edge(-1),
m_nearest_face(-1),
m_prrenderer(NULL)
{
    
}

Sim::~Sim()
{
    
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  General initialization of a simulation
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool Sim::init(const std::string & option_file, bool save_outputs, bool wo_visualization,const std::string env_map_path,const std::string inputdata_dir)
{
    
    // declare and load the options
    Options::addStringOption ("scene", "T1");
    Options::addStringOption ("load-dir", "");
    Options::addDoubleOption ("time-step", 0.01);
    Options::addDoubleOption ("simulation-time", 1.0);
    Options::addBooleanOption("implicit_scheme", false);
    Options::addBooleanOption("RK4-velocity-integration", false);
    Options::addDoubleOption ("smoothing-coef", 0.0);
    Options::addDoubleOption ("damping-coef", 1.0);
    Options::addBooleanOption ("sparse", false);
    
    Options::addBooleanOption ("logging_geometry", false);
    Options::addBooleanOption ("write_geometry", true);
    Options::addBooleanOption ("logging_time", true);
    Options::addBooleanOption ("logging_detailed_time", false);
    
    Options::addBooleanOption("output-png", true);
    Options::addIntegerOption("output-png-every-n-frames", 0);     // 0 means synching with simulation frame rate (equivalent to 1).
    Options::addBooleanOption("output-mesh", false);
    Options::addIntegerOption("output-mesh-every-n-frames", 0);    // 0 means synching with simulation frame rate (equivalent to 1).
    Options::addBooleanOption("output-obj", false);
    Options::addIntegerOption("output-obj-every-n-frames", 0);     // 0 means synching with simulation frame rate (equivalent to 1).
    Options::addBooleanOption("output-ply", false);
    Options::addIntegerOption("output-ply-every-n-frames", 0);     // 0 means synching with simulation frame rate (equivalent to 1).

    Options::addDoubleOption ("remeshing-resolution", 0.1);
    Options::addIntegerOption("remeshing-iterations", 1);
    
    Options::addDoubleOption ("lostopos-collision-epsilon-fraction", 1e-4);       // lostopos collision epsilon (fraction of mean edge length)
    Options::addDoubleOption ("lostopos-merge-proximity-epsilon-fraction", 0.02); // lostopos merge proximity epsilon (fraction of mean edge length)
    Options::addBooleanOption("lostopos-perform-smoothing", false);               // whether or not to perform smoothing
    Options::addDoubleOption ("lostopos-max-volume-change-fraction", 1e-4);       // maximum allowed volume change during a remeshing operation (fraction of mean edge length cubed)
    Options::addDoubleOption ("lostopos-min-triangle-angle", 3.0);                // min triangle angle (in degrees)
    Options::addDoubleOption ("lostopos-max-triangle-angle", 177.0);              // max triangle angle (in degrees)
    Options::addDoubleOption ("lostopos-large-triangle-angle-to-split", 160.0);   // threshold for large angles to be split
    Options::addDoubleOption ("lostopos-min-triangle-area-fraction", 0.02);       // minimum allowed triangle area (fraction of mean edge length squared)
    Options::addBooleanOption("lostopos-t1-transition-enabled", true);            // whether t1 is enabled
    Options::addDoubleOption ("lostopos-t1-pull-apart-distance-fraction", 0.1);   // t1 pull apart distance (fraction of mean edge legnth)
    Options::addBooleanOption("lostopos-smooth-subdivision", false);              // whether to use smooth subdivision during remeshing
    Options::addBooleanOption("lostopos-allow-non-manifold", true);               // whether to allow non-manifold geometry in the mesh
    Options::addBooleanOption("lostopos-allow-topology-changes", true);           // whether to allow topology changes
    
    Options::addIntegerOption("num_subdivision", 0);
    Options::addIntegerOption("mesh-size-n", 2);
    Options::addIntegerOption("mesh-size-m", 2);
    
    Options::addBooleanOption("auto-burst", false);
    Options::addDoubleOption("auto-burst-interval", 10.0);
    Options::addDoubleOption("auto-burst-start", 10.0);
    
    Options::addBooleanOption ("save_mesh", false);
    Options::addBooleanOption ("with_gravity", false);
    Options::addDoubleOption ("gravity_scale", 0.02);//0.02
    
    Options::addBooleanOption ("add_velocity", false);
    Options::addBooleanOption ("add_circular_velocity", false);

    Options::addBooleanOption ("accel", true);
    Options::addDoubleOption ("surface_tension", 1.0);
    Options::addDoubleOption ("gas_pressure", 1.0);

    Options::addStringOption ("sub_scene", "");
    
    Options::addBooleanOption ("two_side", false);
    Options::addBooleanOption ("advect", false);
    
    Options::addBooleanOption ("pulling", false);
    
    Options::addStringOption ("direction", "horizontal");
    
#ifdef FEOS
    Options::addStringOption("vertex_advection", "standard");
    Options::addStringOption("triangle_advection", "standard");
    Options::addBooleanOption ("with_evaporation", false);
    Options::addDoubleOption ("evaporation_per_second", 10.0);
    Options::addDoubleOption ("flow_tension", 0.1);
    Options::addDoubleOption ("vorticity_confinement", 0.0);
    Options::addDoubleOption ("ext_force", 1.0);
    Options::addDoubleOption ("vertical_gravity_scale", 0.02);
    Options::addDoubleOption ("tangential_gravity_scale", 1e-6);//1.0
    Options::addDoubleOption ("thickness_aware_deformation", 0.0);
    Options::addDoubleOption ("deform_standard_thickness", -1.0);
    Options::addDoubleOption ("min_th_nm", 200);
    Options::addDoubleOption ("max_th_nm", 1200);
    Options::addDoubleOption ("burst_ratio", 0.2);
    Options::addDoubleOption ("absorb_liquid", 1.0);//0.02
    Options::addBooleanOption("preserve_local_volumes",false);
    Options::addDoubleOption ("bottom-level", -1000);
    Options::addDoubleOption ("vertical-velocity-coef", 0.1);

#endif //FEOS
    
    Options::parseOptionFile(option_file, m_verbose);
    
    // select the scene
    m_scene = Options::strValue("scene");
    
    if (save_outputs)
    {
        std::stringstream output_dir_ss;
        output_dir_ss << "output_" << ::time(NULL);
        m_output_directory = output_dir_ss.str();
        
#ifndef WIN32
        mkdir( m_output_directory.c_str(), 0755 );
#else
        _mkdir( m_output_directory.c_str() );
#endif // !WIN32
        
        std::cout << "Outputing to directory: " << m_output_directory << std::endl;
    }
    
    if (m_scene == "load")
    {
        m_load_directory = Options::strValue("load-dir");
        assert(m_load_directory != "");
        
        HGF * vs = new HGF(std::vector<LosTopos::Vec3d>(), std::vector<LosTopos::Vec3st>(), std::vector<LosTopos::Vec2i>());
        MeshIO::load(*vs, m_load_directory + "/mesh000000.rec");
        
        std::vector<LosTopos::Vec3d> vertices = vs->m_st->pm_positions;
        std::vector<LosTopos::Vec3st> faces = vs->m_st->m_mesh.m_tris;
        std::vector<LosTopos::Vec2i> face_labels = vs->m_st->m_mesh.m_triangle_labels;
        
        hgf = new HGF(vertices, faces, face_labels);
        MeshIO::load(*hgf, m_load_directory + "/mesh000000.rec");
    } else
    {
        std::vector<LosTopos::Vec3d> vertices;
        std::vector<LosTopos::Vec3st> faces;
        std::vector<LosTopos::Vec2i> face_labels;
        std::vector<size_t> constrained_vertices;
        std::vector<Vec3d> constrained_positions;
        
        if (m_scene == "sphere")                    hgf = Scenes::sceneSphere              (this, vertices, faces, face_labels, constrained_vertices, constrained_positions);
        else if (m_scene == "tet")                  hgf = Scenes::sceneTet                 (this, vertices, faces, face_labels, constrained_vertices, constrained_positions);
        else if (m_scene == "cube")                 hgf = Scenes::sceneCube                (this, vertices, faces, face_labels, constrained_vertices, constrained_positions);
        else if (m_scene == "sheet")                hgf = Scenes::sceneSheet               (this, vertices, faces, face_labels, constrained_vertices, constrained_positions);
        else if (m_scene == "barrel")               hgf = Scenes::sceneBarrel              (this, vertices, faces, face_labels, constrained_vertices, constrained_positions);
        else if (m_scene == "doublebubble")         hgf = Scenes::sceneDoubleBubble        (this, vertices, faces, face_labels, constrained_vertices, constrained_positions);
        else if (m_scene == "twobubbles")           hgf = Scenes::sceneTwoBubbles          (this, vertices, faces, face_labels, constrained_vertices, constrained_positions);
        else if (m_scene == "triplejunction")       hgf = Scenes::sceneTripleJunction      (this, vertices, faces, face_labels, constrained_vertices, constrained_positions);
        else if (m_scene == "foaminit")             hgf = Scenes::sceneFoamInit            (this, vertices, faces, face_labels, constrained_vertices, constrained_positions);
        
        else if (m_scene == "quadjunction")         hgf = Scenes::sceneQuadJunction        (this, vertices, faces, face_labels, constrained_vertices, constrained_positions);
        else if (m_scene == "constrainedsphere")    hgf = Scenes::sceneConstrainedSphere   (this, vertices, faces, face_labels, constrained_vertices, constrained_positions);
        else if (m_scene == "bubblewand")           hgf = Scenes::sceneBubbleWand          (this, vertices, faces, face_labels, constrained_vertices, constrained_positions);
        else if (m_scene == "tworingspinching")     hgf = Scenes::sceneTwoRingsPinching    (this, vertices, faces, face_labels, constrained_vertices, constrained_positions);
        else if (m_scene == "peanutbubble")         hgf = Scenes::scenePeanutBubble        (this, vertices, faces, face_labels, constrained_vertices, constrained_positions);
        else if (m_scene == "straw")                hgf = Scenes::sceneStraw               (this, vertices, faces, face_labels, constrained_vertices, constrained_positions);
        else if (m_scene == "carousel")             hgf = Scenes::sceneCarousel            (this, vertices, faces, face_labels, constrained_vertices, constrained_positions);
        else if (m_scene == "octahedron")           hgf = Scenes::sceneOctahedron          (this, vertices, faces, face_labels, constrained_vertices, constrained_positions);
#ifdef FEOS
        else if (m_scene == "cylinder")             hgf = Scenes::sceneCylinder          (this, vertices, faces, face_labels, constrained_vertices, constrained_positions);
#endif //FEOS
        else if (m_scene == "inputmesh" or m_scene=="inputrec"){
            
            hgf = Scenes::sceneInputData(this, vertices, faces, face_labels, constrained_vertices, constrained_positions,inputdata_dir);
        }
        
        else if (m_scene == "cubicframe")           hgf = Scenes::sceneCubicFrame (this, vertices, faces, face_labels, constrained_vertices, constrained_positions);
        else if (m_scene == "lattice")              hgf = Scenes::sceneLattice (this, vertices, faces, face_labels, constrained_vertices, constrained_positions);
        else if (m_scene == "sphereOverFilm")       hgf = Scenes::sceneSphereOverFilm (this, vertices, faces, face_labels, constrained_vertices, constrained_positions);
        else if (m_scene == "squareFilm")           hgf = Scenes::sceneSquareFilm (this, vertices, faces, face_labels, constrained_vertices, constrained_positions);
        else if (m_scene == "twistedFilm90")        hgf = Scenes::sceneTwistedFilm90 (this, vertices, faces, face_labels, constrained_vertices, constrained_positions);
        else if (m_scene == "brakke")               hgf = Scenes::sceneBrakke (this, vertices, faces, face_labels, constrained_vertices, constrained_positions,inputdata_dir);
        
        if(Options::boolValue("pulling")){
            Scenes::setPulling( hgf);
        }
        
        std::cout << "initial_nv = " << vertices.size() << ", initial_nf = " << faces.size() << std::endl;
        
    }
    
    if(hgf->add_velocity){
        Scenes::set_initial_bubble_velocities(this, hgf);
    }
    
    m_time = 0;
    m_dt = Options::doubleValue("time-step");
    
    m_finished = false;
    
    if(hgf->logging_geometry){
        hgf->geometric_information(m_time);
    }
    
    if(hgf->save_mesh){
        // Sanity check. Terminates the simulation when directory outputmesh does not exist.

        std::string outputmesh_dir="outputmesh";
        struct stat info;
        if((stat( outputmesh_dir.c_str(), &info ) != 0) or not (info.st_mode & S_IFDIR )){
            std::cout<<"Cannot find or access \"outputmesh\" directory.\nPlease make sure you have accessible \"outputmesh\" directory in the right location\nThe proper location is your working directory (On linux and mac, mostly same location as the binary. On windows please check working directory in your Visual Studio)";
            exit(-1);
        }
        outputmesh_dir+="/";
        hgf->write_mesh(outputmesh_dir);
        //hgf->writeMesh_FaceLabel_constrainedVertices();
    }
    
    // output the initial frame
    if (m_output_directory != "" && Options::boolValue("output-mesh"))
        MeshIO::save(*hgf, m_output_directory + "/mesh000000.rec");
    
    // PR rendering
    if (!wo_visualization)
        m_prrenderer = new PRRenderer(hgf,env_map_path);
    
    return true;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  Time stepping
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void Sim::step()
{
    assert(m_scene != "unspecified");
    assert(hgf);
    assert(hgf->surfTrack());
    if (m_verbose)
        std::cout << "Time stepping: t = " << m_time << ", dt = " << m_dt << std::endl;
    
    if(hgf->do_stepScenes){
        
        // scene-specific time stepping
        if (m_scene == "sphere")                    Scenes::stepSphere              (m_dt, this, hgf);
        else if (m_scene == "tet")                  Scenes::stepTet                 (m_dt, this, hgf);
        else if (m_scene == "cube")                 Scenes::stepCube                (m_dt, this, hgf);
        else if (m_scene == "sheet")                Scenes::stepSheet               (m_dt, this, hgf);
        else if (m_scene == "barrel")               Scenes::stepBarrel              (m_dt, this, hgf);
        else if (m_scene == "doublebubble")         Scenes::stepDoubleBubble        (m_dt, this, hgf);
        else if (m_scene == "twobubbles")           Scenes::stepTwoBubbles          (m_dt, this, hgf);
        else if (m_scene == "triplejunction")       Scenes::stepTripleJunction      (m_dt, this, hgf);
        else if (m_scene == "foaminit"){             Scenes::stepFoamInit            (m_dt, this, hgf);
            
        }
        
        else if (m_scene == "quadjunction")         Scenes::stepQuadJunction        (m_dt, this, hgf);
        else if (m_scene == "constrainedsphere")    Scenes::stepConstrainedSphere   (m_dt, this, hgf);
        else if (m_scene == "bubblewand")           Scenes::stepBubbleWand          (m_dt, this, hgf);
        else if (m_scene == "pullingfoam")          Scenes::pullBubbles         (m_dt, this, hgf);
        else if (m_scene == "peanutbubble")         Scenes::stepPeanutBubble        (m_dt, this, hgf);
        else if (m_scene == "straw")                Scenes::stepStraw               (m_dt, this, hgf);
        else if (m_scene == "carousel")             Scenes::stepCarousel            (m_dt, this, hgf);
        else if (m_scene == "octahedron")           Scenes::stepOctahedron          (m_dt, this, hgf);
        else if (m_scene == "tworingspinching")     Scenes::stepTwoRingsPinching    (m_dt, this, hgf);
        
        else if (m_scene == "original_mesh")        Scenes::stepTwoRingsPinching    (m_dt, this, hgf);
        else if (m_scene == "cubicframe")           Scenes::stepCubicFrame    (m_dt, this, hgf);
        else if (m_scene == "sphereOverFilm")       Scenes::stepSphereOverFilm    (m_dt, this, hgf);
        else if (m_scene == "squareFilm")           Scenes::stepSquareFilm    (m_dt, this, hgf);
        else if (m_scene == "twistedFilm90")        Scenes::stepTwistedFilm90    (m_dt, this, hgf);
        
        if (Options::boolValue("pulling")){
            Scenes::pullBubbles(m_dt, this, hgf);
        }
        
    }
    
    if(Options::doubleValue("bottom-level")!=-1000.0){
        Scenes::wall_collision(m_dt, this,hgf);
    }
    
    if(m_scene=="brakke"){
        Scenes::stepBrakke    (m_dt, this, hgf);
        
    }
    
    static bool advect=Options::boolValue("advect");
    if(advect){
        Scenes::move_vertices(m_dt, this,hgf);
    }
    
    if(hgf->add_velocity){
        Scenes:: add_velocities_to_bubbles(m_dt,this, hgf);
    }
    
    if(Options::boolValue("add_circular_velocity")){
        Scenes::add_circular_velocity(m_dt, this, hgf);
    }
    
    // Burst a bubble only when there is more than one bubble.
    if(hgf->auto_burst and hgf->num_currentBubbles >1 ){
        switch (hgf->burstType) {
            case HGF::Random:
                Scenes::burstBubblesRandomly(m_dt, this, hgf);
                break;
            case HGF::Thickness:
                Scenes::burstOneThinBubble(m_dt, this, hgf);
                break;
        }
    }
    
    Scenes::moveLeftOrRight(m_dt, this, hgf);
    Scenes::volume_change(m_dt, this, hgf);
    
    double dt;
    static double computational_time=0.0;
    static double average_computational_time=0.0;
    static int count=0;
    const int count_offset=1;
    
    
    if(hgf->logging_time ){
        
        if(count>=count_offset){
            // general time stepping
            static std::chrono::system_clock::time_point  start, end;
            
            start = std::chrono::system_clock::now();
            
            dt = hgf->step(m_dt);
            
            end = std::chrono::system_clock::now();
            double elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count(); //Convert time to ms.
            
            computational_time+=elapsed;
            average_computational_time=computational_time/(count-count_offset+1);
            
        }else{
            dt = hgf->step(m_dt);
        }
        
    }
    else{
        
        dt= hgf->step(m_dt);
    }
    
    // advance time
    m_frameid++;
    
    m_time+=dt;
    
    if (m_time >= Options::doubleValue("simulation-time"))
        m_finished = true;
    
    if(hgf->logging_geometry){
        hgf->geometric_information(m_time);
    }
    if(hgf->logging_time){
        
        if(count>=count_offset){
            std::cout<<"average_time:"<<average_computational_time<<std::endl;
            
            if(hgf->logging_detailed_time){
                std::cout<<"dAdx_time:"<<hgf->time_dAdx<<",vc_time:"<<hgf->time_volume_correction<<",los_topos_time:"<<hgf->time_los_topos<<",other_time:"<<computational_time-(hgf->time_dAdx+hgf->time_volume_correction+hgf->time_los_topos)<<std::endl;
            }
        }
        
        count++;
        
    }
    
    if(hgf->save_mesh){
        
        const double record_interval=0.001;
        const double eps=1e-6;
        static double passed_time=0;
        passed_time+=dt;
        
        if(passed_time>record_interval-eps){
            passed_time=0;
            
            hgf->write_mesh("outputmesh/");
            //hgf->writeMesh_FaceLabel_constrainedVertices();
        }
        
    }
    
}
#ifdef HAVE_PNG
void Sim::outputOneImage(){
    
    std::cout<<"Saving image to image.png...";
    int w, h;
    w = glutGet(GLUT_WINDOW_WIDTH);
    h = glutGet(GLUT_WINDOW_HEIGHT);
    
    YImage img;
    img.resize(w, h);
    glReadPixels(0, 0, w, h, GL_RGBA, GL_UNSIGNED_BYTE, (unsigned char *)(img.data()));
    img.flip();
    img.save("image.png");
    std::cout<<"done.\n";
}
#endif

void Sim::stepOutput(bool wo_visualization)
{
    
    if (m_output_directory != "")
    {
        
        int frameid = (int)(time() / dt() + 0.5);
        
#ifdef HAVE_PNG 
        int pngfd = Options::intValue("output-png-every-n-frames");
        if ((pngfd == 0 || frameid % pngfd == 0) && Options::boolValue("output-png") && !wo_visualization)
        {
            static int image_count=0;
            
            std::stringstream png_ss;
            
            png_ss << m_output_directory << "/frame" << std::setfill('0') << std::setw(6) << ++image_count << ".png";
            
            int w, h;
            w = glutGet(GLUT_WINDOW_WIDTH);
            h = glutGet(GLUT_WINDOW_HEIGHT);
            
            YImage img;
            img.resize(w, h);
            glReadPixels(0, 0, w, h, GL_RGBA, GL_UNSIGNED_BYTE, (unsigned char *)(img.data()));
            img.flip();
            img.save(png_ss.str().c_str());
        }
#endif //HAVE_PNG
        
        int meshfd = Options::intValue("output-mesh-every-n-frames");
        if ((meshfd == 0 || frameid % meshfd == 0) && Options::boolValue("output-mesh"))
        {
            static int rec_count=0;
            
            std::stringstream mesh_ss;
            mesh_ss << m_output_directory << "/state" << std::setfill('0') << std::setw(6) << rec_count++ << ".rec";
            MeshIO::save(*hgf, mesh_ss.str());
        }
        
        int objfd = Options::intValue("output-obj-every-n-frames");
        if ((objfd == 0 || frameid % objfd == 0) && Options::boolValue("output-obj"))
        {
            static int obj_count=0;
            std::stringstream obj_ss;
            obj_ss << m_output_directory << "/mesh" << std::setfill('0') << std::setw(6) << obj_count++ << ".obj";
            const bool with_imaginary_vertices=false;
            MeshIO::saveOBJ(*hgf, obj_ss.str(),with_imaginary_vertices,hgf->saveMeshWith_normal,hgf->saveMeshChange_y_z);
        }
        
        int plyfd = Options::intValue("output-ply-every-n-frames");
        if ((plyfd == 0 || frameid % plyfd == 0) && Options::boolValue("output-ply"))
        {
            static int ply_count=0;
            std::stringstream ply_ss;
            ply_ss << m_output_directory << "/mesh" << std::setfill('0') << std::setw(6) << ply_count++ << ".ply";
            const bool with_imaginary_vertices=false;
            MeshIO::savePLY(*hgf, ply_ss.str(),with_imaginary_vertices,hgf->saveMeshChange_y_z);
        }
        
    }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  Loading saved simulation
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool Sim::load(int inc)
{
    int current_frame = (int)((m_time + m_dt * 0.5) / m_dt);
    int next_frame = current_frame + inc;
    if (next_frame < 0) next_frame = 0;
    
    std::stringstream ss;
    ss << m_load_directory << "/mesh" << std::setfill('0') << std::setw(6) << next_frame << ".rec";
    if (!MeshIO::load(*hgf, ss.str()))
    {
        std::cout << "Loading frame " << ss.str() << " unsuccessful." << std::endl;
        return false;
    }
    
    m_time = m_dt * next_frame;
    
    return true;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  Rendering
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
namespace
{
    bool is_edge_nonmanifold(const LosTopos::SurfTrack & st, size_t e)
    {
        return st.m_mesh.m_edge_to_triangle_map[e].size() != 2;
    }
    
    bool is_vertex_nonmanifold(const LosTopos::SurfTrack & st, size_t v)
    {
        for (size_t i = 0; i < st.m_mesh.m_vertex_to_edge_map[v].size(); i++)
            if (is_edge_nonmanifold(st, st.m_mesh.m_vertex_to_edge_map[v][i]))
                return true;
        return false;
    }
    
    bool is_face_next_to_nonmanifold_vertices(const LosTopos::SurfTrack & st, size_t f)
    {
        return is_vertex_nonmanifold(st, st.m_mesh.m_tris[f][0]) || is_vertex_nonmanifold(st, st.m_mesh.m_tris[f][1]) || is_vertex_nonmanifold(st, st.m_mesh.m_tris[f][2]);
    }
    
    bool is_edge_next_to_nonmanifold_vertices(const LosTopos::SurfTrack & st, size_t e)
    {
        return is_vertex_nonmanifold(st, st.m_mesh.m_edges[e][0]) || is_vertex_nonmanifold(st, st.m_mesh.m_edges[e][1]);
    }
    
    bool is_vertex_next_to_nonmanifold_vertices(const LosTopos::SurfTrack & st, size_t v)
    {
        for (size_t i = 0; i < st.m_mesh.m_vertex_to_edge_map[v].size(); i++)
            if (is_edge_next_to_nonmanifold_vertices(st, st.m_mesh.m_vertex_to_edge_map[v][i]))
                return true;
        return false;
    }
}

void Sim::render(RenderMode rm, const Vec2d & mousepos, int selection_mask)
{
    if (rm == RM_PR)
    {
        if(m_prrenderer->env_map_path==""){
            return;
        }
        m_prrenderer->render();
        return;
    }
    
    // find the primitive being picked by cursor
    Mat4d MV;
    Mat4d PJ;
    {
        float mv[16];
        glGetFloatv(GL_MODELVIEW_MATRIX, mv);
        float pj[16];
        glGetFloatv(GL_PROJECTION_MATRIX, pj);
        MV << mv[0], mv[4], mv[8], mv[12], mv[1], mv[5], mv[9], mv[13], mv[2], mv[6], mv[10], mv[14], mv[3], mv[7], mv[11], mv[15];
        PJ << pj[0], pj[4], pj[8], pj[12], pj[1], pj[5], pj[9], pj[13], pj[2], pj[6], pj[10], pj[14], pj[3], pj[7], pj[11], pj[15];
        
    }
    
    //Display camera information for Mitsuba renderer.
    if(camera_information){
        camera_information=false;
        
        Eigen::Vector3d lookat_up(MV(1,0),MV(1,1),MV(1,2));
        Eigen::Vector3d lookat_target(-MV(2,0),-MV(2,1),-MV(2,2));
        Eigen::Vector3d lookat_origin;
        
        double back_scale=2.;
        for(int i=0;i<3;++i){
            lookat_origin[i]=-back_scale*(MV(0,i)*MV(0,3)+MV(1,i)*MV(1,3)+MV(2,i)*MV(2,3));
        }
        
        std::cout<<"<lookat ";
        std::cout<<"target=\"";
        for(int i=0;i<3;++i){
            std::cout<<lookat_target[i]<<" ";
        }
        std::cout<<"\" ";
        
        std::cout<<"origin=\"";
        for(int i=0;i<3;++i){
            std::cout<<lookat_origin[i]<<" ";
        }
        std::cout<<"\" ";
        
        std::cout<<"up=\"";
        for(int i=0;i<3;++i){
            std::cout<<lookat_up[i]<<" ";
        }
        std::cout<<"\"";
        std::cout<<"/>";
        std::cout<<std::endl;
    }
    
    Mat4d MVP = PJ * MV;
    
    double mind = -1;
    m_nearest_vertex = -1;
    m_nearest_edge = -1;
    m_nearest_face = -1;
    if (selection_mask & SM_VERTEX)
    {
        for (size_t i = 0; i < hgf->mesh().nv(); i++)
        {
            Vec3d pos = hgf->pos(i);
            Vec4d scrpos_h = MVP * Vec4d(pos.x(), pos.y(), pos.z(), 1.0);
            Vec2d scrpos = Vec2d(scrpos_h.x(), scrpos_h.y()) / scrpos_h.w();
            
            double distance = (scrpos - mousepos).norm();
            if (distance < mind || mind < 0)
            {
                mind = distance;
                m_nearest_vertex = i;
            }
        }
    }
    
    if (selection_mask & SM_EDGE)
    {
        for (size_t i = 0; i < hgf->mesh().ne(); i++)
        {
            Vec3d v0 = hgf->pos(hgf->mesh().m_edges[i][0]);
            Vec3d v1 = hgf->pos(hgf->mesh().m_edges[i][1]);
            
            Vec4d scrv0_h = MVP * Vec4d(v0.x(), v0.y(), v0.z(), 1.0);
            Vec2d scrv0 = Vec2d(scrv0_h.x(), scrv0_h.y()) / scrv0_h.w();
            Vec4d scrv1_h = MVP * Vec4d(v1.x(), v1.y(), v1.z(), 1.0);
            Vec2d scrv1 = Vec2d(scrv1_h.x(), scrv1_h.y()) / scrv1_h.w();
            
            double distance = (mousepos - (scrv0 + scrv1) / 2).norm();
            //            double distance = (mousepos - scrv0 - (mousepos - scrv0).dot(scrv1 - scrv0) * (scrv1 - scrv0) / (scrv1 - scrv0).squaredNorm()).norm();
            if (distance < mind || mind < 0)
            {
                mind = distance;
                m_nearest_vertex = -1;
                m_nearest_edge = i;
            }
        }
    }
    
    if (selection_mask & SM_FACE)
    {
        for (size_t i = 0; i < hgf->mesh().nt(); i++)
        {
            const LosTopos::Vec3st & t = hgf->mesh().get_triangle(i);
            Vec3d v0 = hgf->pos(t[0]);
            Vec3d v1 = hgf->pos(t[1]);
            Vec3d v2 = hgf->pos(t[2]);
            
            Vec4d scrv0_h = MVP * Vec4d(v0.x(), v0.y(), v0.z(), 1.0);
            Vec2d scrv0 = Vec2d(scrv0_h.x(), scrv0_h.y()) / scrv0_h.w();
            Vec4d scrv1_h = MVP * Vec4d(v1.x(), v1.y(), v1.z(), 1.0);
            Vec2d scrv1 = Vec2d(scrv1_h.x(), scrv1_h.y()) / scrv1_h.w();
            Vec4d scrv2_h = MVP * Vec4d(v2.x(), v2.y(), v2.z(), 1.0);
            Vec2d scrv2 = Vec2d(scrv2_h.x(), scrv2_h.y()) / scrv2_h.w();
            
            double distance = (mousepos - (scrv0 + scrv1 + scrv2) / 3).norm();
            if (distance < mind || mind < 0)
            {
                mind = distance;
                m_nearest_vertex = -1;
                m_nearest_edge = -1;
                m_nearest_face = i;
            }
        }
    }
    
    assert(mind >= 0);
    assert(m_nearest_vertex >= 0 || m_nearest_edge >= 0 || m_nearest_face >= 0);
    
    bool truncate = false;
    
    glEnable(GL_DEPTH_TEST);
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    glShadeModel(GL_SMOOTH);
    glEnable(GL_POLYGON_OFFSET_FILL);
    glPolygonOffset(1.0, 1.0);
    
    if (rm == RM_TRANSPARENT || rm == RM_NONMANIFOLD)
    {
        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        glDepthMask(GL_FALSE);
    } else
    {
        glEnable(GL_LIGHTING);
        glEnable(GL_LIGHT0);
        
        GLfloat mat_diffuse[] = { 1.0, 1.0, 1.0, 1.0 };
        glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
        GLfloat mat_specular[] = { 1.0, 1.0, 1.0, 1.0 };
        glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
        GLfloat mat_shininess[] = { 50.0 };
        glMaterialfv(GL_FRONT, GL_SHININESS, mat_shininess);
        
        GLfloat light_ambient[] = { 0.3, 0.3, 0.3, 1.0 };
        glLightfv(GL_LIGHT0, GL_AMBIENT, light_ambient);
        GLfloat light_direction[] = { 1.0, 1.0, 1.0, 0.0 };
        glLightfv(GL_LIGHT0, GL_SPOT_DIRECTION, light_direction);
    }
    
    // pre-compute vertex normals (area-weighted face normals), for RM_OPAQUE_SMOOTH_SHADED mode
    std::vector<Vec3d> vn(hgf->mesh().nv(), Vec3d(0, 0, 0));
    for (size_t i = 0; i < hgf->mesh().nt(); i++)
    {
        LosTopos::Vec3st t = hgf->surfTrack()->m_mesh.get_triangle(i);
        Vec3d x0 = hgf->pos(t[0]);
        Vec3d x1 = hgf->pos(t[1]);
        Vec3d x2 = hgf->pos(t[2]);
        
        Vec3d nt = (x1 - x0).cross(x2 - x0);
        if (hgf->surfTrack()->m_mesh.get_triangle_label(i)[0] < hgf->surfTrack()->m_mesh.get_triangle_label(i)[1])
            nt = -nt;
        
        vn[t[0]] += nt;
        vn[t[1]] += nt;
        vn[t[2]] += nt;
    }
    for (size_t i = 0; i < hgf->mesh().nv(); i++)
        vn[i].normalize();
    
    if (false)
    {
        glLineWidth(5);
        glBegin(GL_LINES);
        for (size_t i = 0; i < hgf->mesh().nt(); i++)
        {
            LosTopos::Vec3st t = hgf->surfTrack()->m_mesh.get_triangle(i);
            Vec3d x0 = hgf->pos(t[0]);
            Vec3d x1 = hgf->pos(t[1]);
            Vec3d x2 = hgf->pos(t[2]);
            Vec3d c = (x0 + x1 + x2) / 3;
            Vec3d n = (x1 - x0).cross(x2 - x0).normalized();
            Vec3d eout = c + n * 0.03;
            Vec3d ein = c - n * 0.03;
            
            LosTopos::Vec2i l = hgf->mesh().get_triangle_label(i);
            
            if (l[1] == 0)      glColor3d(1, 0, 0);
            else if (l[1] == 1) glColor3d(0, 1, 0);
            else if (l[1] == 2) glColor3d(0, 0, 1);
            else if (l[1] == 3) glColor3d(0.8, 0.8, 0);
            else if (l[1] == 4) glColor3d(0.8, 0, 0.8);
            else if (l[1] == 5) glColor3d(0, 0.8, 0.8);
            else if (l[1] == 6) glColor3d(0.4, 0.4, 1);
            else if (l[1] == 7) glColor3d(0.4, 1, 0.4);
            else if (l[1] == 8) glColor3d(1, 0.4, 0.4);
            else                glColor3d(0, 0, 0);
            glVertex3d(c[0], c[1], c[2]);
            glVertex3d(eout[0], eout[1], eout[2]);
            
            if (l[0] == 0)      glColor3d(1, 0, 0);
            else if (l[0] == 1) glColor3d(0, 1, 0);
            else if (l[0] == 2) glColor3d(0, 0, 1);
            else if (l[0] == 3) glColor3d(0.8, 0.8, 0);
            else if (l[0] == 4) glColor3d(0.8, 0, 0.8);
            else if (l[0] == 5) glColor3d(0, 0.8, 0.8);
            else if (l[0] == 6) glColor3d(0.4, 0.4, 1);
            else if (l[0] == 7) glColor3d(0.4, 1, 0.4);
            else if (l[0] == 8) glColor3d(1, 0.4, 0.4);
            else                glColor3d(0, 0, 0);
            glVertex3d(c[0], c[1], c[2]);
            glVertex3d(ein[0], ein[1], ein[2]);
        }
        glEnd();
        glLineWidth(1);
        
    }
    
    // Render triangles
    glBegin(GL_TRIANGLES);
    
    for (size_t i = 0; i < hgf->surfTrack()->m_mesh.nt(); i++)
    {
        if (rm == RM_NONMANIFOLD)
        {
            if (!is_face_next_to_nonmanifold_vertices(*hgf->surfTrack(), i))
                continue;
        }
        
        if (hgf->m_st->triangle_is_all_solid(i))
            continue;
        
        LosTopos::Vec3st t = hgf->surfTrack()->m_mesh.get_triangle(i);
        LosTopos::Vec3d x0 = hgf->surfTrack()->pm_positions[t[0]];
        LosTopos::Vec3d x1 = hgf->surfTrack()->pm_positions[t[1]];
        LosTopos::Vec3d x2 = hgf->surfTrack()->pm_positions[t[2]];
        LosTopos::Vec3d c = (x0 + x1 + x2) / 3;
        
        if (truncate)
            if (c[0] > 0.5 || c[0] < -0.5)
                continue;
        
        double shrink = (rm == RM_TRANSPARENT ? 0.05 : 0);
        x0 += (c - x0) * shrink;
        x1 += (c - x1) * shrink;
        x2 += (c - x2) * shrink;
        Vec3d n0, n1, n2;
        
        if (rm == RM_OPAQUE_FLAT_SHADED)
        {
            n0 = vc(cross(x1 - x0, x2 - x0));
            if (hgf->surfTrack()->m_mesh.get_triangle_label(i)[0] < hgf->surfTrack()->m_mesh.get_triangle_label(i)[1])
                n0 = -n0;
            n0.normalize();
            n1 = n0;
            n2 = n0;
        } else
        {
            n0 = vn[t[0]];
            n1 = vn[t[1]];
            n2 = vn[t[2]];
        }
        
        if (m_nearest_face == i)
            glColor4d(0.4, 0.5, 0.6, 0.5);
        else
            // Triangle color
            glColor4d(0.7, 0.8, 0.9, 0.0);
        
        glNormal3d(n0[0], n0[1], n0[2]);    glVertex3d(x0[0], x0[1], x0[2]);
        glNormal3d(n1[0], n1[1], n1[2]);    glVertex3d(x1[0], x1[1], x1[2]);
        glNormal3d(n2[0], n2[1], n2[2]);    glVertex3d(x2[0], x2[1], x2[2]);
    }
    glEnd();
    
    if (rm == RM_TRANSPARENT || rm == RM_NONMANIFOLD)
    {
        glDisable(GL_BLEND);
        glDepthMask(GL_TRUE);
    } else
    {
        glDisable(GL_LIGHTING);
        glDisable(GL_LIGHT0);
    }
    
    // render edges
    if (rm != RM_OPAQUE_SMOOTH_SHADED)
    {
        glLineWidth(1);
        glColor3d(0.3, 0.3, 0.3);
        glBegin(GL_LINES);
        for (size_t i = 0; i < hgf->mesh().ne(); i++)
        {
            if(hgf->unrendered_edges[i]){
                continue;
            }
            if (rm == RM_NONMANIFOLD)
            {
                //                if (!is_edge_next_to_nonmanifold_vertices(*hgf->surfTrack(), i))
                if (!is_edge_nonmanifold(*hgf->surfTrack(), i))
                    continue;
            }
            
            Vec3d x0 = hgf->pos(hgf->mesh().m_edges[i][0]);
            Vec3d x1 = hgf->pos(hgf->mesh().m_edges[i][1]);
            
            if (truncate)
                if (x0.x() + x1.x() > 1 || x0.x() + x1.x() < -1)
                    continue;
            
            glVertex3d(x0[0], x0[1], x0[2]);
            glVertex3d(x1[0], x1[1], x1[2]);
        }
        glEnd();
    }
    
#ifdef FEOS
    
        size_t num_v=hgf->surfTrack()->m_mesh.nv();
        size_t num_e=hgf->edge_thickness.size();
        size_t num_t=hgf->triangle_centers.size();
    
// Give colors to triangles
    if( hgf->facePaint != HGF::FacePaint::NONE){
            
            double alpha_tri_vel=0.3;
            if (rm == RM_TRANSPARENT || rm == RM_NONMANIFOLD)
            {
                // this flag determines transparency
                if (hgf->triangles_transparent) {
                    glEnable( GL_BLEND );
                    glBlendFunc( GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA );
                    glDepthMask( GL_FALSE );
                }
                
                alpha_tri_vel = 1.0; 0.4;
                
                
            }else{
                alpha_tri_vel=1.0;
            }
            
            glBegin(GL_TRIANGLES);

			Eigen::MatrixXd Colors;
            hgf->colorForRendering(Colors,hgf->draw_th_power);
            
            std::vector<Vec3d> vn(hgf->mesh().nv(), Vec3d(0, 0, 0));
            const size_t num_rt=hgf->real_triangles.size();
            for (size_t rti = 0; rti <num_rt; rti++)
            {
                const size_t i = hgf->real_triangles[rti];
                
                if(hgf->unrendered_triangles[i]){continue;}

                LosTopos::Vec3st t = hgf->surfTrack()->m_mesh.get_triangle(i);
                
                for (size_t iv = 0; iv < 3; ++iv) {
                    size_t vert = t[ iv ];
                    LosTopos::Vec3d x = hgf->surfTrack()->pm_positions[ vert ];
                    Eigen::Vector4d c;
                    c<<Colors( vert, 0 ), Colors( vert, 1 ),  Colors( vert, 2 ), 0.2;

                    glColor4d(c[ 0 ], c[ 1 ], c[ 2 ], c[3] );
                    glVertex3d( x[ 0 ], x[ 1 ], x[ 2 ] );
                    
                }
                
            }
            glEnd();
            
            glPointSize(2);
            glBegin(GL_LINES);
            
            //velocity on surface
            //3d velocity
            const bool display_3d_vertex_velocity=0;
            if(display_3d_vertex_velocity){
                for (size_t i = 0; i <num_v; i++)
                {
                    if (rm == RM_NONMANIFOLD)
                    {
                        //                if (!is_vertex_next_to_nonmanifold_vertices(*hgf->surfTrack(), i))
                        if (!is_vertex_nonmanifold(*hgf->surfTrack(), i))
                            continue;
                    }
                    
                    const LosTopos::Vec3d & x = hgf->surfTrack()->pm_positions[i];
                    
                    if (truncate)
                        if (x[0] > 0.5 || x[0] < -0.5)
                            continue;
                    glColor3f(1.f, 0, 0);
                    glVertex3d(x[0], x[1], x[2]);
                    
                    const auto& velocity=hgf->vel(i);
                    const double vel_scale=1.0;
                    
                    auto dest= x +vc(vel_scale*velocity);
                    glColor3f(0, 1.f, 0);
                    glVertex3d(dest[0], dest[1], dest[2]);
                    
                    
                    
                }
            }
            
            for(size_t e=0;e< num_e;++e){
                
                size_t v0=hgf->mesh().m_edges[e][0];
                size_t v1=hgf->mesh().m_edges[e][1];
                
                Eigen::Vector3d x=get_hgf()->edge_centers[e];
                
                const bool display_edge_3D_velocity=0;
                if(display_edge_3D_velocity){
                    glColor3f(1.f, 0, 0);
                    glVertex3d(x[0], x[1], x[2]);
                    const auto& velocity0=hgf->vel(v0);
                    const auto& velocity1=hgf->vel(v1);
                    const auto& velocity=(velocity0+velocity1)/2.0;
                    
                    const double edge_3D_vel_scale=1.0;
                    
                    auto dest = x +edge_3D_vel_scale*velocity;
                    glColor3f(0, 1.f, 0);
                    glVertex3d(dest[0], dest[1], dest[2]);
                }
                
                const bool display_edge_normal=0;
                if(display_edge_normal){
                    glColor3f(1.f, 0, 0);
                    glVertex3d(x[0], x[1], x[2]);
                    
                    
                    const auto & normal=hgf->edge_t_normals[e];
                    const double normal_scale=0.1;//.04;//0.02;
                    auto dest=x +normal_scale*normal;
                    glColor3f(0, 1.f, 0);
                    glVertex3d(dest[0], dest[1], dest[2]);
                }
                
            }
            
            // Draw flow velocities
            if (hgf->draw_vertex_velocity) {
                for (size_t i = 0; i <num_v; i++)
                {
                    //if (rm == RM_NONMANIFOLD)
                    //{
                    //    //                if (!is_vertex_next_to_nonmanifold_vertices(*hgf->surfTrack(), i))
                    //    if (!is_vertex_nonmanifold( *hgf->surfTrack(), i ))
                    //        continue;
                    //}
                    
                    const LosTopos::Vec3d & x = hgf->surfTrack()->pm_positions[ i ];
                    
                    if (truncate)
                        if (x[ 0 ] > 0.5 || x[ 0 ] < -0.5)
                            continue;
                    
                    glColor3f( 1.f, 0, 0 );
                    glVertex3d( x[ 0 ], x[ 1 ], x[ 2 ] );
                    
                    const auto& velocity = hgf->fv( i );
                    
                    auto dest = x + vc( hgf->draw_vector_scale*velocity );
                    glColor3f( 0, 1.f, 0 );
                    glVertex3d( dest[ 0 ], dest[ 1 ], dest[ 2 ] );
                    
                    
                }
                
            }
            
            bool display_pressure_grad_vertex=0;
            if (display_pressure_grad_vertex and hgf->GradP_vert.size()>0) {
                for (size_t i = 0; i <num_v; i++)
                {
                    const LosTopos::Vec3d & x = hgf->surfTrack()->pm_positions[ i ];
                    
                    glColor3f( 0.f, 0, 0 );
                    
                    glVertex3d( x[ 0 ], x[ 1 ], x[ 2 ] );
                    
                    const auto& gradP = hgf->GradP_vert.row(i);
                    
                    auto dest = x + vc( hgf->draw_vector_scale*gradP );
                    glColor3f( 0, 0.f, 1.f );
                    glVertex3d( dest[ 0 ], dest[ 1 ], dest[ 2 ] );
                }
            }
            
            for(size_t t=0;t<num_t;++t){
                
                if(hgf->draw_triangle_velocity){
                    glColor3f(1.f, 0, 0);
                    const auto&tri_center=get_hgf()->triangle_centers[t];
                    glVertex3d(tri_center[0], tri_center[1], tri_center[2]);
                    
                    const auto & tri_velocity=hgf->triangle_velocities[t];
                    auto dest=tri_center +hgf->draw_vector_scale*tri_velocity;
                    glColor3f(0, 1.f, 0);
                    glVertex3d(dest[0], dest[1], dest[2]);
                }
                
            }
            
// For debugging backtrace in advection
//#define tracedraw
#ifdef tracedraw
            if(hgf->intersectionOnce){
                glColor3f(0,0.5f, 1.f);
                auto &trace=hgf->trace;
                for (size_t i = 0; i < trace.size()-1; i++)
                {
                    glVertex3d(trace[i][0], trace[i][1], trace[i][2]);
                    glVertex3d(trace[i+1][0], trace[i+1][1], trace[i+1][2]);
                    
                }
            }
#endif tracedraw
            
            glEnd();
            glPointSize(1);
#endif //FEOS
            
            if (rm == RM_TRANSPARENT || rm == RM_NONMANIFOLD)
            {
                if (hgf->triangles_transparent) {
                    glDisable(GL_BLEND);
                    glDepthMask(GL_TRUE);
                }
                
            }
            
        }

// For debugging backtrace in advection
//#define debug_projected_adj
#ifdef debug_projected_adj
    
    glPointSize(2);
    glColor3f(1, 0, 0);
    glBegin(GL_POINTS);
    
    for(auto& adj_vAndProj:hgf->proj_points_around_0th_ver){
        size_t adj_ind=adj_vAndProj.first;
        
        glVertex3d(adj_vAndProj.second[0], adj_vAndProj.second[1], adj_vAndProj.second[2]);
    }
    
    static bool once=false;
    if(!once){
        for(auto& adj_vAndProj:hgf->proj_points_around_0th_ver){
            size_t adj_ind=adj_vAndProj.first;
            
            once=true;
        }
        
    }
    
    glEnd();
#endif //debug_projected_adj
    
    // render vertices, with mean curvature coloring
    // actually doing nothing as we do not keep mean curvature.
    if (rm != RM_OPAQUE_SMOOTH_SHADED and hgf->draw_auxiliaries)
    {
        glPointSize(2);
        glColor3f(0, 0, 0);
        glBegin(GL_POINTS);
        for (size_t i = 0; i < hgf->surfTrack()->m_mesh.nv(); i++)
        {
            if (rm == RM_NONMANIFOLD)
            {
                //                if (!is_vertex_next_to_nonmanifold_vertices(*hgf->surfTrack(), i))
                if (!is_vertex_nonmanifold(*hgf->surfTrack(), i))
                    continue;
            }
            
            const LosTopos::Vec3d & x = hgf->surfTrack()->pm_positions[i];
            
            if (truncate)
                if (x[0] > 0.5 || x[0] < -0.5)
                    continue;
            
            //glVertex3d(x[0], x[1], x[2]);
        }
        glEnd();
        glPointSize(1);
    }
    
    // render constrained vertices
    if (rm != RM_OPAQUE_SMOOTH_SHADED)
    {
        glPointSize(8);
        glColor3f(0, 0.7, 1);
        glBegin(GL_POINTS);
        for (size_t i = 0; i < hgf->constrainedPositions().size(); i++)
        {
            Vec3d x = hgf->constrainedPositions()[i];
            glVertex3d(x[0], x[1], x[2]);
        }

        glEnd();
        glPointSize(10);
        glColor3f(1, 0.7, 0);
        glBegin(GL_POINTS);
        for (size_t i = 0; i < hgf->mesh().nv(); i++)
        {
            Vec3d x = hgf->pos(i);
            if (hgf->m_st->vertex_is_any_solid(i))
                glVertex3d(x[0], x[1], x[2]);
        }
        glEnd();
        glPointSize(1);
    }
    
    
    
    if (hgf->draw_auxiliaries and m_nearest_vertex >= 0)
    {
        glPointSize(12);
        glColor3f(0, 0, 0);
        glBegin(GL_POINTS);
        Vec3d x = hgf->pos(m_nearest_vertex);
        glVertex3d(x[0], x[1], x[2]);
        glEnd();
    }
    
    if (hgf->draw_auxiliaries and m_nearest_edge >= 0)
    {
        glColor3d(0, 0, 0);
        glLineWidth(3);
        glBegin(GL_LINES);
        Vec3d x0 = hgf->pos(hgf->mesh().m_edges[m_nearest_edge][0]);
        Vec3d x1 = hgf->pos(hgf->mesh().m_edges[m_nearest_edge][1]);
        glVertex3d(x0[0], x0[1], x0[2]);
        glVertex3d(x1[0], x1[1], x1[2]);
        glEnd();
        glLineWidth(1);
    }
    
}

void Sim::showPrimitiveInfo()
{
    if (m_nearest_vertex >= 0)
    {
        std::cout << "Vertex of Interest: " << m_nearest_vertex << " (" << hgf->pos(m_nearest_vertex).transpose() << ")" << std::endl;
        
        std::cout << "  incident edges:"; for (size_t i = 0; i < hgf->mesh().m_vertex_to_edge_map    [m_nearest_vertex].size(); i++) std::cout << " " << hgf->mesh().m_vertex_to_edge_map    [m_nearest_vertex][i]; std::cout << std::endl;
        std::cout << "  incident faces:"; for (size_t i = 0; i < hgf->mesh().m_vertex_to_triangle_map[m_nearest_vertex].size(); i++) std::cout << " " << hgf->mesh().m_vertex_to_triangle_map[m_nearest_vertex][i]; std::cout << std::endl;
        
    }
    
    if (m_nearest_edge >= 0)
    {
        std::cout << "Edge of Interest: " << m_nearest_edge << ": " << hgf->mesh().m_edges[m_nearest_edge][0] << " (" << hgf->pos(hgf->mesh().m_edges[m_nearest_edge][0]).transpose() << ") - " << hgf->mesh().m_edges[m_nearest_edge][1] << " (" << hgf->pos(hgf->mesh().m_edges[m_nearest_edge][1]).transpose() << ") length = " << (hgf->pos(hgf->mesh().m_edges[m_nearest_edge][1]) - hgf->pos(hgf->mesh().m_edges[m_nearest_edge][0])).norm() << std::endl;
        
        std::cout << "  incident faces:"; for (size_t i = 0; i < hgf->mesh().m_edge_to_triangle_map[m_nearest_edge].size(); i++) std::cout << " " << hgf->mesh().m_edge_to_triangle_map[m_nearest_edge][i]; std::cout << std::endl;
    }
    
    if (m_nearest_face >= 0)
    {
        std::cout << "Face of Interest: " << m_nearest_face << ": " << hgf->mesh().m_tris[m_nearest_face][0] << " (" << hgf->pos(hgf->mesh().m_tris[m_nearest_face][0]).transpose() << "), " << hgf->mesh().m_tris[m_nearest_face][1] << " (" << hgf->pos(hgf->mesh().m_tris[m_nearest_face][1]).transpose() << "), " << hgf->mesh().m_tris[m_nearest_face][2] << " (" << hgf->pos(hgf->mesh().m_tris[m_nearest_face][2]).transpose() << ")" << std::endl;
        std::cout << "  labels: " << hgf->mesh().get_triangle_label(m_nearest_face) << std::endl;
        std::cout << "  incident edges:"; for (size_t i = 0; i < 3; i++) std::cout << " " << hgf->mesh().m_triangle_to_edge_map[m_nearest_face][i]; std::cout << std::endl;
    }
}

