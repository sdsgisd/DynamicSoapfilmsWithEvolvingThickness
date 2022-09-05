//
//  eigenheaders.cpp
//  MultiTracker
//
//  Created by Fang Da on 10/27/14.
//
//

#include "eigenheaders.h"

Vec3d vc(const LosTopos::Vec3d & v)
{
    return Vec3d(v[0], v[1], v[2]);
}

LosTopos::Vec3d vc(const Vec3d & v)
{
    return LosTopos::Vec3d(v[0], v[1], v[2]);
}

MatXd vc(const std::vector<LosTopos::Vec3d> & vs)
{
    const size_t num=vs.size();
    MatXd eigenMat;
    eigenMat.resize(num, 3);
    for(size_t i=0;i<num;++i){
        eigenMat.row(i)=vc(vs[i]);
    }
    return eigenMat;
}
