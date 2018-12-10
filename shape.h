//
// Created by mzy on 2018/11/5.
//

#ifndef RENDERING_SHAPE_H
#define RENDERING_SHAPE_H

#include "core.h"
#include "spectrum.h"
#include "BSDF.h"
struct intersection {
    Vec3d nShading, nGeometry, hit, wo;
    spectrum spec;
    double a;
    BSDF *bsdf;
    bool isnull;
    bool operator < (intersection const &x) const {
        return a < x.a;
    }
};


intersection nullIsect();

class shape {
public:
    //virtual bool intersect(ray r, double &a) const = 0;
    virtual Matx<double, 3, 2> minmax() const = 0;
    int intersectPlane(int dim, double pos) {
        auto mm = minmax();
        int flag = 0;
        if (mm(dim, 0) < pos) flag |= 1;
        if (mm(dim, 1) > pos) flag |= 2;
        return flag;
    }
    virtual intersection intersect(ray r) const = 0;
};

class triangle: public shape {

protected:
    Vec3d points[3], nPoints[3], nGeometry;
    bool smoothed;

public:



    triangle(Vec3d points_[3], Vec3d *nPoints_ = nullptr) {
        rep(i, 3) points[i] = points_[i];
        if (nPoints_ == nullptr) smoothed = false;
        else {
            rep(i, 3) nPoints[i] = nPoints_[i];
            smoothed = true;
        }
        nGeometry = (points[1] - points[0]).cross(points[2] - points[0]);
        nGeometry = nGeometry / norm(nGeometry);

    }
    triangle(Vec3d u, Vec3d v, Vec3d w, Vec3d offset = Vec3d(0, 0, 0)) {
        points[0] = u + offset;
        points[1] = v + offset;
        points[2] = w + offset;

        smoothed = false;
        nGeometry = (points[1] - points[0]).cross(points[2] - points[0]);
        nGeometry = nGeometry / norm(nGeometry);

    }
    //bool intersect(ray r, double &a) const;
    Matx<double, 3, 2> minmax() const;
    intersection intersect(ray r) const;
};

class revoBezier: public triangle {

    vector<Vec3d> pointsC;
    double theta;
public:


    revoBezier(vector<Vec3d> pointsC, Vec3d points3[3], double theta): pointsC(pointsC), triangle(points3), theta(theta) {

    }

    //bool intersect(ray r, double &a) const;
    intersection intersect(ray r) const override;
};

revoBezier *subRevoBezier(vector<Vec3d> points, int iTheta, int totTheta, int iT, int totT, int j);

#endif //RENDERING_SHAPE_H
