//  main.cpp
//
//  Created by Fang Da on 2014.
//
//  Modified by Sadashige Ishida in 2017.

#define ADTEST

#include <iostream>
#include <sstream>

#ifdef __APPLE__
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#include <GLUT/glut.h>
#elif __linux__
#include <GL/glew.h>
#include <GL/glut.h>
#include <GL/gl.h>
#include <GL/glu.h>
#elif _WIN32
#include <windows.h>
#include <GL/glew.h>
#include <GL/gl.h>
#include <GL/glu.h>
//#include <GLUT/glut.h>
#include <GL/glut.h>
#endif

#include "Options.h"
#include "Sim.h"
#include "MeshIO.h"

Sim g_sim(false);

struct SimControl
{
    bool wo_visualization;
    
    int win_w;
    int win_h;
    
    bool step;
    bool run;
    bool autoload;
    
    double view_theta;
    double view_alpha;
    double view_dist;
    
    int mouse_x;
    int mouse_y;
    
    bool ldrag;
    int ldrag_start_x;
    int ldrag_start_y;
    bool rdrag;
    int rdrag_start_x;
    int rdrag_start_y;
    Sim::RenderMode render_mode;
    
#ifdef FEOS
    bool mdrag{false};
    int mdrag_start_x;
    int mdrag_start_y;
#endif
    
    int selection_mode;
    
} g_sc;

void renderBitmapString(float x, float y, float z, std::string s)
{
    
    glColor3f(0, 0, 0);
    glRasterPos3f(x, y, z);
    for (size_t i = 0; i < s.size(); i++)
        glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, s[i]);
}

void display()
{
    auto hgf=g_sim.get_hgf();
    // object center and zoom
    Vec3d center(0, 0, 0);
    for (size_t i = 0; i < hgf->mesh().nv(); i++)
        center += hgf->pos(i);
    center /= hgf->mesh().nv();
    Vec3d radius(0, 0, 0);
    for (size_t i = 0; i < hgf->mesh().nv(); i++)
        for (size_t j = 0; j < 3; j++)
            radius[0] = std::max(radius[0], (hgf->pos(i) - center)[0]);
    
    // Limiter for camera view. 
    double min_d = std::max(std::max(radius[0], radius[1]), radius[2]) * 1.2;
    //g_sc.view_dist = std::max(min_d, g_sc.view_dist);
    
    glClearColor(1, 1, 1, 1);
    glClearDepth(1);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    double ar = (double)g_sc.win_w / g_sc.win_h;
    double vfh = 0.01 * 0.4;
    glFrustum(-vfh * ar, vfh * ar, -vfh, vfh, 0.01, 500);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    
    glTranslated(0, 0, -g_sc.view_dist);
    glRotated(-90, 1, 0, 0);
    glRotated(g_sc.view_alpha, 1, 0, 0);
    glRotated(g_sc.view_theta, 0, 0, 1);
        
    if(hgf->draw_auxiliaries){
        glBegin(GL_LINES);
        glColor3d(1, 0, 0);     glVertex3d(0, 0, 0);    glVertex3d(2, 0, 0);
        glColor3d(0, 1, 0);     glVertex3d(0, 0, 0);    glVertex3d(0, 2, 0);
        glColor3d(0, 0, 1);     glVertex3d(0, 0, 0);    glVertex3d(0, 0, 2);
        glEnd();
        
    }
    
    g_sim.render(g_sc.render_mode, Vec2d((double)g_sc.mouse_x / g_sc.win_w * 2 - 1, 1 - (double)g_sc.mouse_y / g_sc.win_h * 2), g_sc.selection_mode);
    
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(0, g_sc.win_w, 0, g_sc.win_h, -1, 1);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    
    if(hgf->draw_auxiliaries){
        std::stringstream ss;
        ss << "T = " << g_sim.time();
        std::string s = ss.str();
        renderBitmapString(20, g_sc.win_h - 20, 0, s);
    }
    glutSwapBuffers();
}

void idle()
{
    if (g_sc.run || g_sc.step)
    {
        g_sc.step = false;
        g_sim.step();
        std::cout << "Finished step: T = " << g_sim.time() << std::endl;
        g_sim.stepOutput(g_sc.wo_visualization);
        if (g_sim.isFinished()){
            g_sim.get_hgf()->writeMesh_FaceLabel_constrainedVertices();
            exit(0);
        }
        
        if (!g_sc.wo_visualization)
            glutPostRedisplay();
    }
    
    if (g_sc.autoload)
    {
        if (!g_sim.load(1))
            exit(0);
        
        g_sim.stepOutput(g_sc.wo_visualization);
        
        if (!g_sc.wo_visualization)
            glutPostRedisplay();
    }
}

void keyboard(unsigned char k, int x, int y)
{
    using namespace Eigen;
    auto hgf=g_sim.get_hgf();
    size_t num_t = hgf->mesh().nt();
    size_t num_e = hgf->mesh().ne();
    size_t num_v = hgf->mesh().nv();
    
    enum Element{
        vertex,
        triangle
    };
    Element element = triangle;
    
    switch(k )
    {
        case 27:
        case 'q':
        case 'Q':
            exit( 0 );
        case ' ':
            g_sc.run ^= 1;
            break;
        case 's':
            g_sc.step = true;
            break;
        case 'D':{
            hgf->do_stepScenes^=1;
            break;
        }
        case'm':
            g_sc.render_mode = ( Sim::RenderMode )( ( (int)g_sc.render_mode + ( k == 'm' ? 1 : -1 ) ) % ( (int)Sim::RM_COUNT ) );
            std::cout << "Render mode: " << (int)g_sc.render_mode << std::endl;
            glutPostRedisplay();
            break;
// Some key settings from SoapFilm3D disabled
//        case 'V':
//            g_sc.selection_mode = ( k == 'v' ? ( g_sc.selection_mode | Sim::SM_VERTEX ) : ( g_sc.selection_mode & ~Sim::SM_VERTEX ) );
//            std::cout << "Mouse cursor selecting" << ( ( g_sc.selection_mode & Sim::SM_VERTEX ) ? " vertices" : "" ) << ( ( g_sc.selection_mode & Sim::SM_EDGE ) ? " edges" : "" ) << ( ( g_sc.selection_mode & Sim::SM_FACE ) ? " faces" : "" ) << "." << std::endl;
//            glutPostRedisplay();
//            break;
//        case 'E':
//            g_sc.selection_mode = ( k == 'e' ? ( g_sc.selection_mode | Sim::SM_EDGE ) : ( g_sc.selection_mode & ~Sim::SM_EDGE ) );
//            std::cout << "Mouse cursor selecting" << ( ( g_sc.selection_mode & Sim::SM_VERTEX ) ? " vertices" : "" ) << ( ( g_sc.selection_mode & Sim::SM_EDGE ) ? " edges" : "" ) << ( ( g_sc.selection_mode & Sim::SM_FACE ) ? " faces" : "" ) << "." << std::endl;
//            glutPostRedisplay();
//            break;
//        case 'F':
//            g_sc.selection_mode = ( k == 'f' ? ( g_sc.selection_mode | Sim::SM_FACE ) : ( g_sc.selection_mode & ~Sim::SM_FACE ) );
//            std::cout << "Mouse cursor selecting" << ( ( g_sc.selection_mode & Sim::SM_VERTEX ) ? " vertices" : "" ) << ( ( g_sc.selection_mode & Sim::SM_EDGE ) ? " edges" : "" ) << ( ( g_sc.selection_mode & Sim::SM_FACE ) ? " faces" : "" ) << "." << std::endl;
//            glutPostRedisplay();
//            break;
//        case 'n':
//        case 'N':
//            g_sim.showPrimitiveInfo();
//            break;
        case 'o':
        case 'O':{
			hgf->geometric_information(g_sim.m_dt);
            const bool write_imaginary_vertices = k == 'o' ? false : true;
            
            hgf->writeMesh_FaceLabel_constrainedVertices( write_imaginary_vertices );
            hgf->write_constrained_mesh( "./constrained_mesh.obj" );
            
#ifdef HAVE_PNG
            g_sim.outputOneImage();
#endif 
            break;
        }
        case '+':
            //Increase the volume of the 0-th bubble.
            hgf->blowing_bubble0 = !g_sim.get_hgf()->blowing_bubble0;
            break;
        case '-':
            //Decrease the volume of the 0-th bubble.
            hgf->absorbing_bubble0 = !g_sim.get_hgf()->absorbing_bubble0;
            break;
        case 'b':{ // Manually burst a bubble
            
#ifdef FEOS
            int burstedRegion=-1;
            if(hgf->burstType==HGF::Random){
                //Burst a randomly chosen bubble.cone
                burstedRegion=Scenes::burstOneBubble(g_sim.m_dt, &g_sim, hgf);
            }
            if(hgf->burstType==HGF::Thickness){
                //Burst a randomly chosen bubble.cone
                burstedRegion=Scenes::burstThinnestBubble(g_sim.m_dt, &g_sim, hgf);
                std::cout<<"region "<<burstedRegion<<" bursted.\n";
            }
            if(burstedRegion!=-1){
                hgf->test_update_mesh_via_LosTopos();
                hgf->init_FES();
            }else{
                std::cout<<"no region bursted.\n";
            }
#else //FEOS
            Scenes::burstOneBubble(g_sim.m_dt, &g_sim, hgf);
#endif //FEOS
            break;
        }
        case 'R':
            //Move constrained vertices to right.
            hgf->move_right = !g_sim.get_hgf()->move_right;
            hgf->move_left = false;
            break;
        case 'L':
            //Move constrained vertices to left.
            hgf->move_left = !g_sim.get_hgf()->move_left;
            hgf->move_right = false;
            break;
        case 'l':
            g_sim.camera_information = true;
            break;
#pragma mark Keys For FEOS
#ifdef FEOS
        case 'i':{ // Change simulation type
            
            hgf->simulationType =static_cast<HGF::SimulationType>( ( hgf->simulationType + 1 ) % hgf->simulationType_names.size() );
            
            hgf->init_FES_preserveVelocities();
            std::cout << "Simulation Type " <<hgf->simulationType<<". "<<hgf->simulationType_names[hgf->simulationType] <<std::endl;
            break;
        }
        
        case 'h':{ // Display current thickness
            Eigen::VectorXd ThTri;
            hgf->scalarVert2triNaive(hgf->thvv(), ThTri);
            const double currentLowestTh=ThTri.minCoeff();
            const double currentLowestThRatio=(currentLowestTh-hgf->min_th_nm)/(hgf->max_th_nm-hgf->min_th_nm);
            std::cout<<"Lowest th = "<<currentLowestTh<<": Ratio="<<currentLowestThRatio<<"\n";
            
            break;
        }
            
        case 'g':{ // Toggle gravity
            hgf->with_gravity^=1;
            std::cout<<"with_gravity:"<<hgf->with_gravity<<"\n";

            break;
        }
        case 'e': // Toggle evaporation
            hgf->with_evaporation^=1;
            std::cout<<"with_evaporation:"<<hgf->with_evaporation<<"\n";
            break;
        case 'B':{ // Set the ratio of thickness triggering bubble burst.
            
            hgf->auto_burst=true;
            double burst_ratio=(hgf->burst_threshold-hgf->min_th_nm)/(hgf->max_th_nm-hgf->min_th_nm);
            std::cout<<"Current burst_ratio="<<burst_ratio<<"\n";
            std::cout<<"Set a new burst_ratio.\nburst_ratio=";
            std::cin>>burst_ratio;
            hgf->burst_threshold=hgf->min_th_nm+burst_ratio*(hgf->max_th_nm-hgf->min_th_nm);
            
            break;
        }
            
        case '|':{ // reset camera position
            g_sc.view_alpha=0.0;
            g_sc.view_theta = 180;
            break;
        }
        case 'A':{
            hgf->draw_auxiliaries^=1;
            break;
        }
            
        case '.':{
            if(hgf->facePaint==HGF::FacePaint::THICKNESS){
                hgf->facePaint=HGF::FacePaint::NONE;
                std::cout<<"HGF::FacePaint=NONE\n";
            }
            else if(hgf->facePaint==HGF::FacePaint::NONE){
                hgf->facePaint=HGF::FacePaint::THICKNESS;
                std::cout<<"HGF::FacePaint=THICKNESS\n";
            }
            break;
        }
            
        case '\\':{
            hgf->triangles_transparent^=1;
            break;
        }

        case 'u':{ // Make a bubble region transparent
            ++hgf->unrendered_region;
            if(hgf->unrendered_region==hgf->num_region){
                hgf->unrendered_region=-1;
            }
            hgf->setUnrenderedTrianglesAndEdges();
            break;
        }
        case 'n':{ // normalize rendered bubble thickness for very subtle thickness distributions
            hgf->normalize_color^=1;
            std::cout<<"color is normalized?"<<hgf->normalize_color<<"\n";
            break;
        }
        case 'v':
            hgf->draw_vertex_velocity ^= 1;
            break;
        case 't':
            hgf->draw_triangle_velocity ^= 1;
            break;
        case '9':
            hgf->draw_th_power *= 0.9;
            std::cout << "draw_th_power=" << hgf->draw_th_power << std::endl;
            break;
        case '0':
            hgf->draw_th_power *= 1.1;
            std::cout << "draw_th_power=" << hgf->draw_th_power << std::endl;
            break;
        case '{':
            hgf->draw_vector_scale *= 5.0;
            break;
        case '}':
            hgf->draw_vector_scale *= 0.2;
            break;

        case '/':{ // Change equation between our model and nabier stokes (Stable Fluids)
            hgf->equation= static_cast<HGF::Equation>( ( hgf->equation + 1 ) % hgf->equation_names.size() );
            
            glutSetWindowTitle( hgf->equation_names[hgf->equation].c_str());
            hgf->init_FES();

            break;
        }
        case 'r':
            //reset velocities and scalar values
            
            hgf->triangle_velocities = std::vector<Eigen::Vector3d>( num_t, Eigen::Vector3d::Zero() );
            hgf->edge_tnormal_velocities = std::vector<double>( num_e, 0.0 );
            
            hgf->triangle_thickness = std::vector<double>( num_t, 0.0);
            hgf->edge_thickness = std::vector<double>( num_e, 0.0 );
            
            hgf->DivV_vert.setZero();
            hgf->DivV_tri.setZero();
            for( int vi = 0; vi < num_v; ++vi ) {
                
                hgf->fv( vi ) <<0.0, 0.0, 0.0 ;
                hgf->thv( vi ) =  hgf->default_th;
                hgf->uthv( vi ) = 0.0;

                hgf->vel(vi)<<0,0,0;
                
            }
            for(int ti=0;ti<num_t;++ti){
                hgf->ft(ti)<<0.0,0.0,0.0;
            }
           
            hgf->init_FES();

            break;

        case 'k':{ // Give thickness to a few vertices
            int partition = hgf->mesh().nv() / 50;
            static int ipk = 0;
            
            int num_target_verts = 1;
            bool withNeighbors = false;
            
            for( int i = 0; i < num_target_verts; ++i ){
                if( ++ipk == partition ){ ipk = 0; }
                int ivert = num_v * ipk / partition;
                
                double sign = ipk % 2 == 0 ? -1 : 1;
                const float thickness_coef=1.0;
                
                hgf->thv( ivert ) += sign *thickness_coef* (hgf->max_th_nm-hgf->min_th_nm)/2.0;
                if( withNeighbors ){
                    std::vector<size_t> adj_verts;
                    hgf->mesh().get_adjacent_vertices( ivert, adj_verts );
                    for( size_t adj_vert : adj_verts ){
                        auto& adj_th=hgf->thv( adj_vert ) ;
                        adj_th+= 0.3*sign *thickness_coef* (hgf->max_th_nm-hgf->min_th_nm)/2.0;
                        
                    }
                }
                
                for(size_t vi=0;vi<num_v;++vi){
                    auto& th=hgf->thv( vi);
                    th=HGF::clamp(th,hgf->min_th_nm,hgf->max_th_nm);
                }
                
            }
            
            hgf->initialLiquidVolume=hgf->computeTotalLiquidVolume();
            break;
        }
        case 'K':{ // Single velocity to a vertex/triangle/edge

            
            int partition = hgf->mesh().nv() / 55; 50; 10;
            static int ipk = 0;
            
            int num_target_elements= 1;
            bool withNeighbors = false;
                        
            for( int i = 0; i < num_target_elements; ++i ){
                
                switch( element ) {
                    case vertex:
                    {
                        size_t num_v = hgf->mesh().nv();
                        if( ++ipk == partition ){ ipk = 0; }
						int ivert = num_v * ipk / partition;
                        auto& fv = hgf->fv( ivert );
                        
                        fv += hgf->ext_force*Vec3d( 200, 4, 4 );
                        
                        const auto vn = vc( hgf->surfTrack()->get_vertex_normal( ivert ) );
                        fv -= fv.dot( vn ) / vn.squaredNorm()*vn;
                        
                        if( withNeighbors ){
                            std::vector<size_t> adj_verts;
                            hgf->mesh().get_adjacent_vertices( ivert, adj_verts );
                            for( size_t adj_vert : adj_verts ){
                                auto& fv_adj = hgf->fv( adj_vert );
                                fv_adj += Vec3d( 60, 4, 4 );
                                const auto vn_adj = vc( hgf->surfTrack()->get_vertex_normal( adj_vert ) );
                                fv_adj -= fv_adj.dot( vn_adj ) / vn.squaredNorm()*fv_adj;
                            }
                        }
                        hgf->velocityVertex2Triangle();
                        
                        break;
                    }
                    case triangle:{
                        // triangle
                        size_t itri = num_t * ipk++ / std::max(partition,1);

                        const double coef=hgf->ext_force;
                        hgf->triangle_velocities[ itri ] +=
                        Vec3d( coef*100, 0.01, 0.01 );

                        const auto tn = hgf->triangle_normals[itri];
                        hgf->triangle_velocities[ itri ] -= hgf->triangle_velocities[ itri ].dot( tn ) / tn.squaredNorm()*tn;
                        hgf->velocityTriangle2Vertex();
                        break;
                    }
                    default:
                        break;
                }
            }
            hgf->setft(hgf->getTriVelocity());
            break;
        }

        case 'p':{ // peform projection
            
                    hgf->project_tri( g_sim.m_dt, hgf->equation);
                    hgf->velocityTriangle2Vertex();

            hgf->setft(hgf->getTriVelocity());
        }

        case 'a':{ // peform advection of thickness
            
            const int num_iteration = 10;
            for( int ite = 0; ite < num_iteration; ++ite ){
                
                    VectorXd Th=hgf->thvv();
                    hgf->advect_vertex( g_sim.m_dt, Th );
                    hgf->setThv(Th);
                    hgf->velocityTriangle2Vertex(true);
                
                g_sim.m_time += g_sim.m_dt;
                hgf->writeVelocityOfPoint( g_sim.m_time, true );
                
            }
            break;
        }
            
        // Add gradually varying height field.
        case 'G':{
            
            Scenes::addThicknessForGradientTest( g_sim.m_dt, &g_sim, hgf );
            
            enum LaplacianMethod{
                // As Direct's result is the half of Direct, we multiply Divgrad by 1/2.
                // Then the both yield the same result.
                DivGrad,
                Direct
            };
            
            LaplacianMethod laplacianMethod=Direct;
            
            MatrixXd TempLT;
            
            switch(laplacianMethod){
                case DivGrad:{
                    Eigen::VectorXd T(num_v);
                    for(size_t vi=0;vi<num_v;++vi){
                        T[vi]=hgf->thv(vi);
                    }

                    hgf->vert2triGrad(T,hgf->GradT_tri);
                    
                    // Compute div grad T
                    hgf->computeWeightedDivergenceOnVertices(hgf->GradT_tri, hgf->LT_vert);

                    break;
                    
                }
                    
                case Direct:{

                    hgf->th_laplacian_filter(.0); //for
                    break;
                }
                    
            }
            hgf->initialLiquidVolume=hgf->computeTotalLiquidVolume();
            break;
            
        }

        case 'X':{ // Add ring-like velocity field.
            
            for(size_t ti=0; ti<num_t; ti++){
                const double y_upper = 0.8;
                const double y_lower = 0.5;
                
                const auto tc = hgf->triangle_centers[ti];
                if((y_lower <= tc.y()) && (tc.y() <= y_upper)) {
                    
                    const Vector3d tn = hgf->triangle_normals[ti];
                    const Vector3d y_axis =Vector3d(0.0, 1.0, 0.0);
                    
                    const double coef=5.0;
                    hgf->triangle_velocities[ ti ] += coef * tn.cross(y_axis);
                    
                    hgf->triangle_velocities[ ti ] -= hgf->triangle_velocities[ ti ].dot( tn ) / tn.squaredNorm()*tn;
                }
            }
            hgf->setft(hgf->getTriVelocity());
            hgf->velocityTriangle2Vertex();
            break;
        }
        case 'P':{ // Add thickness randomly by Perlin noise
            int dimension;
            std::cout << "int dimension = ";
            std::cin >> dimension;
            
            double frequency;
            std::cout << "double frequency = ";
            std::cin >> frequency;
            
            int octaves;
            std::cout << "int32 octaves    = ";
            std::cin >> octaves;
            
            std::uint32_t seed;
            std::cout << "uint32 seed      = ";
            std::cin >> seed;
            Scenes::addThicknessByPerlinNoise( g_sim.m_dt, &g_sim, hgf, dimension, frequency, seed, octaves );
            break;
        }

        case 'C':{ // Add flow velocity by curl noise.

            double frequency;
            std::cout << "double frequency = ";
            std::cin >> frequency;
            
            int octaves;
            std::cout << "int32 octaves    = ";
            std::cin >> octaves;
            
            std::uint32_t seed;
            std::cout << "uint32 seed      = ";
            std::cin >> seed;
            Scenes::curlNoise( g_sim.m_dt, &g_sim, hgf, frequency, seed, octaves );
            hgf->setft(hgf->getTriVelocity());
            break;
        }
        case '*':{ // initialize thickness and velocity by noise
            Scenes::addThicknessByPerlinNoise( g_sim.m_dt, &g_sim, hgf, 1, 20, 1, 1 );
            Scenes::addThicknessByPerlinNoise( g_sim.m_dt, &g_sim, hgf, 1, 29, 41, 0 );
            Scenes::curlNoise( g_sim.m_dt, &g_sim, hgf, 21, 1245, 0 );
            hgf->setft(hgf->getTriVelocity());
            break;
        }

        case ':':{
            Scenes::addVerticalVelocity(g_sim.m_dt,&g_sim,hgf);
            break;
        }
            
#endif //FEOS
    }
    
}

void mouse(int b, int s, int x, int y)
{
    if (b == GLUT_LEFT_BUTTON && s == GLUT_DOWN)
    {
        g_sc.ldrag = true;
        g_sc.ldrag_start_x = x;
        g_sc.ldrag_start_y = y;
    } else if (b == GLUT_LEFT_BUTTON && s == GLUT_UP)
    {
        g_sc.ldrag = false;
    } else if (b == GLUT_RIGHT_BUTTON && s == GLUT_DOWN)
    {
        g_sc.rdrag = true;
        g_sc.rdrag_start_x = x;
        g_sc.rdrag_start_y = y;
    } else if (b == GLUT_RIGHT_BUTTON && s == GLUT_UP)
    {
        g_sc.rdrag = false;
    }
#ifdef FEOS
    else if (b == GLUT_MIDDLE_BUTTON && s == GLUT_DOWN)
    {
        g_sc.mdrag = true;
        g_sc.mdrag_start_x = x;
        g_sc.mdrag_start_y = y;
    }else if (b == GLUT_MIDDLE_BUTTON && s == GLUT_UP)
    {
        g_sc.mdrag = false;
    }
#endif
    
    glutPostRedisplay();
}

void motion(int x, int y)
{
    auto hgf=g_sim.get_hgf();
    if (g_sc.ldrag)
    {
        g_sc.view_theta += (x - g_sc.ldrag_start_x) * 1.0;
        g_sc.view_alpha += (y - g_sc.ldrag_start_y) * 1.0;
        
        g_sc.ldrag_start_x = x;
        g_sc.ldrag_start_y = y;
    }
    if (g_sc.rdrag)
    {
        g_sc.view_dist *= pow(2.0, (y - g_sc.rdrag_start_y) * 0.01);
        
        g_sc.rdrag_start_x = x;
        g_sc.rdrag_start_y = y;
    }
    
    g_sc.mouse_x = x;
    g_sc.mouse_y = y;
    
#ifdef FEOS
    if(g_sc.mdrag){
        
        // next account for direction using delta_alpha, (delta_beta)
        double delta_x=(x - g_sc.mdrag_start_x) * 1e-8;
        double delta_y=(y - g_sc.mdrag_start_y) * 1e-8;
        
        Eigen::Quaterniond theta_rotate( Eigen::AngleAxisd(g_sc.view_theta, Eigen::Vector3d::UnitZ()));
        Eigen::Quaterniond alpha_rotate( Eigen::AngleAxisd(g_sc.view_alpha, Eigen::Vector3d::UnitX()));
        Eigen::Quaterniond vertical_rotate( Eigen::AngleAxisd(-M_PI/4.0, Eigen::Vector3d::UnitX()));
        
        Eigen::Vector3d delta(delta_x,0,delta_y);
        delta=theta_rotate*alpha_rotate*vertical_rotate*delta;
        //delta=vertical_rotate*alpha_rotate*theta_rotate*delta;

        bool vertex_base=0;
        if(vertex_base){
            int vertex=g_sim.m_nearest_vertex;
            if(vertex==-1 ){
                return;
            }
            Eigen::Vector3d vert_normal(0,0,0);
            
            vert_normal=vc(hgf->surfTrack()->get_vertex_normal(vertex));
            
            delta-=delta.dot(vert_normal)*vert_normal/        vert_normal.squaredNorm();
            
            hgf->fv( vertex )+=delta;
            
        }
        
        bool tri_base=0;
        if(tri_base){
            int tri=g_sim.m_nearest_face;
            if(tri==-1 ){
                return;
            }
            
            hgf->triangle_velocities[tri]+=delta;
            
            
        }
        
        bool toNearestVertex=1;
        if(toNearestVertex){
            // Find nearest vertex
            Eigen::Vector3d screen_center(0,0,-g_sc.view_dist);
            
            Eigen::Quaterniond minus_vertical_rotate( Eigen::AngleAxisd(M_PI/4.0, Eigen::Vector3d::UnitX()));
            Eigen::Quaterniond minus_alpha_rotate( Eigen::AngleAxisd(-g_sc.view_alpha, Eigen::Vector3d::UnitX()));
            Eigen::Quaterniond minus_theta_rotate( Eigen::AngleAxisd(-g_sc.view_theta, Eigen::Vector3d::UnitZ()));
            
            screen_center=minus_vertical_rotate*minus_alpha_rotate* minus_theta_rotate*screen_center;
            size_t nv=hgf->mesh().nv();
            double mindist=1e8;
            int nearest_v=-1;
            for(int vi=0;vi<nv;++vi){
                double dist=(screen_center-hgf->pos(vi) ).squaredNorm();
                if (dist<mindist){
                    mindist=dist;
                    nearest_v=vi;
                }
            }
            if(nearest_v==-1 ){
                return;
            }
            Eigen::Vector3d vert_normal(0,0,0);
            
            vert_normal=vc(hgf->surfTrack()->get_vertex_normal(nearest_v));
            
            delta-=delta.dot(vert_normal)*vert_normal/        vert_normal.squaredNorm();
            
            hgf->fv( nearest_v )+=delta;
        }
        
        
        
    }
    
#endif 
    
    glutPostRedisplay();
}

void passiveMotion(int x, int y)
{
    g_sc.mouse_x = x;
    g_sc.mouse_y = y;
    
    glutPostRedisplay();
}

void reshape(int w, int h)
{
    g_sc.win_w = w;
    g_sc.win_h = h;
    
    glViewport(0, 0, w, h);
    
    glutPostRedisplay();
}

void parse_arguments(int argc, char **argv,bool *output_ptr,bool*wo_visualization_ptr,std::string*env_map_path_ptr,std::string*inputdata_dir){
    
    using std::cout;
    using std::endl;
    using std::string;
    using std::stringstream;
    
    for(int i=1;i<argc;++i){
        
        stringstream stream(argv[i]);
        string attribute_name;
        string attribute_value;
        
        getline(stream,attribute_name,'=');
        getline(stream,attribute_value,'=');
        
        //        if(attribute_value.empty()){
        //            continue;
        //        }
        //        assert(!attribute_value.empty());
        
        if(attribute_name=="output"){
            *output_ptr=1;
        }
        else if(attribute_name=="wo_visualization"){
            *wo_visualization_ptr=1;
        }
        else if(attribute_name=="env_map_path"){
            *env_map_path_ptr= attribute_value;
        }
        else if(attribute_name=="inputdata_dir"){
            *inputdata_dir= attribute_value;
        }
        
    }
    
}

int main(int argc, char * argv[])
{
    
    std::cout << "[Usage] specify a configuration file as the command line arguments. e.g. doublebubble.txt" << std::endl;
    std::cout<< "[Other options]"<<std::endl;
    std::cout << "output: for outputting images and meshes." << std::endl;
    std::cout << "wo_visualization: Simulating without visualization." << std::endl;
    std::cout << "env_map_path=\"ENV_MAP_PATH\": Specify the directory of the environment map." << std::endl;
    std::cout << "inputdata_dir=\"INPUT_DATA_DIR\": Specify the directory of the input data such as obj or rec." << std::endl;
    std::cout<<std::endl;
    
    if(argc==1){exit(1);}
    
    // simulation setup
    g_sc.run = false;
    g_sc.step = false;
    g_sc.autoload = false;
    
    g_sc.win_w = 640;//1280;
    g_sc.win_h = 720;
    
    g_sc.view_theta = 180;
    g_sc.view_alpha = 0;// Using 180 makes the scene upside-down.
    g_sc.view_dist = 4;
    
    g_sc.mouse_x = 0;
    g_sc.mouse_y = 0;
    
    g_sc.ldrag = false;
    g_sc.ldrag_start_x = 0;
    g_sc.ldrag_start_y = 0;
    g_sc.rdrag = false;
    g_sc.rdrag_start_x = 0;
    g_sc.rdrag_start_y = 0;
    
    g_sc.render_mode = Sim::RM_TRANSPARENT;
    
    g_sc.selection_mode = Sim::SM_VERTEX | Sim::SM_EDGE | Sim::SM_FACE;
    
    if (!g_sc.wo_visualization)
    {
        // glut setup
        glutInit(&argc, argv);
#ifdef _WIN32
        glutInitWindowPosition( 2*g_sc.win_w, g_sc.win_h );
#else
        glutInitWindowPosition( 0,0 );
#endif
        glutInitWindowSize(g_sc.win_w, g_sc.win_h);
        glutCreateWindow("Soap Film Equation");
        
        glutKeyboardFunc(keyboard);
        glutMouseFunc(mouse);
        glutMotionFunc(motion);
        glutPassiveMotionFunc(passiveMotion);
        glutIdleFunc(idle);
        glutDisplayFunc(display);
        glutReshapeFunc(reshape);
    }
    
    bool output=false;
    std::string env_map_path="";
    std::string inputdata_dir="";
    
    parse_arguments(argc, argv, &output,&(g_sc.wo_visualization),&env_map_path,&inputdata_dir);
    
    bool success = g_sim.init(argv[1], output, g_sc.wo_visualization,env_map_path,inputdata_dir);
    if (!success)
        return 1;
    
    std::cout << "Initialization complete. Starting the simulation..." << std::endl;
    
    // main loop
    if (g_sc.wo_visualization)
    {
        g_sc.run = true;
        while (true)
            idle();
    } else
    {
        glutMainLoop();
    }
    
}

