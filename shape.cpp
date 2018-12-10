//
// Created by mzy on 2018/11/5.
//

#include "shape.h"



intersection nullIsect() {
    intersection isect;
    isect.isnull = true;
    return isect;
}

intersection triangle::intersect(ray r) const {
    intersection ret;
    Matx<double, 3, 3> A;
    rep(i, 3) rep(j, 3) A(i, j) = points[j][i] - r.o[i];
    if (abs(determinant(A)) < eps) return nullIsect();
    auto coef = A.inv() * r.dir;
    double a = 0;
    rep(i, 3) {
        if (coef(i) < eps) return nullIsect();
        a += coef(i);
    }
    ret.hit = r.o + r.dir / a;
    ret.nGeometry = nGeometry;
    ret.nShading = Vec3d();
    if (smoothed) {
        rep(i, 3) ret.nShading += coef(i) * nPoints[i];
        ret.nShading /= a;
    }
    else ret.nShading = ret.nGeometry;
    if (ret.nShading.dot(r.dir) > 0) ret.nShading *= -1;


    ret.nShading /= norm(ret.nShading);
    ret.a = a;
    ret.wo = -r.dir / norm(r.dir);

    ret.isnull = false;
    return ret;
}

/*
bool triangle::intersect(ray r, double &a) const {

    Matx<double, 3, 3> A;
    rep(i, 3) rep(j, 3) A(i, j) = points[j][i] - r.o[i];
    auto toInv = A;
    if (abs(determinant(toInv)) < eps) return false;
    auto coef = toInv.inv() * r.dir;
    a = 0;
    rep(i, 3) {
        if (coef(i) < eps) return false;
        a += coef(i);
    }
    return true;
}
*/
Matx<double, 3, 2> triangle::minmax() const {
    double minx[3], maxx[3];
    Matx<double, 3, 2> ret;
    Matx<double, 3, 3> tmp;
    rep(i, 3) rep(j, 3) tmp(i, j) = points[j](i);
    rep(i, 3) minMaxIdx(tmp.row(i), &ret(i, 0), &ret(i, 1));
    return ret;
}

Vec3d deCast(vector<Vec3d> points, double t) {
    auto n = points.size();
    for (int i = n - 1; i >= 0; i--)
        rep(j, i)
            points[j] = (1 - t) * points[j] + t * points[j + 1];
    return points[0];
}

revoBezier *subRevoBezier(vector<Vec3d> points, int iTheta, int totTheta, int iT, int totT, int j) {
    double theta1 = 2 * acos(-1)  * iTheta / totTheta, theta2 = 2 * acos(-1) * (iTheta + 1) / totTheta;
    double t1 = (double)iT / totT, t2 = (iT + 1.0) / totT;
    Vec3d point0 = deCast(points, t1), point1 = deCast(points, t2), points3[3];
    auto n = points.size() - 1;
    int index = 0;
    points3[0] = Vec3d(point0[0] * cos(theta1), point0[0] * sin(theta1), point0[2]);
    if (j == 0) points3[++index] = Vec3d(point0[0] * cos(theta2), point0[0] * sin(theta2), point0[2]);
    if (j == 1) points3[++index] = Vec3d(point1[0] * cos(theta1), point1[0] * sin(theta1), point1[2]);
    points3[++index] = Vec3d(point1[0] * cos(theta2), point1[0] * sin(theta2), point1[2]);
    if (j == 0) swap(points3[1], points3[2]);
    return new revoBezier(points, points3, (theta1 + theta2) / 2);
}


intersection revoBezier::intersect(ray r) const {
    //return triangle::intersect(r);
    auto ret = triangle::intersect(r); //new intersection;
    double theta, t, tt;
    if (!ret.isnull) {
        theta = this->theta;
        t = 0.5;
        tt = 1 / ret.a;
    }
    else {
        theta = this->theta;
        t = 0.5;
        tt = norm((points[0] + points[1] + points[2]) / 3 - r.o) / norm(r.dir);
    }

    auto n = pointsC.size() - 1;
    Matx<double, 3, 3> J;
    Vec3d res;
    int max_iter = 20;

    for (;;) {

        if (--max_iter < 0) break;
        auto cur = deCast(pointsC, t);
        res[0] = cos(theta) * cur[0];
        res[1] = sin(theta) * cur[0];
        res[2] = cur[2];
        rep(i, 3) res[i] -= r.o[i] + tt * r.dir[i];

        vector<Vec3d> tmp;
        for (int i = 1; i <= n; i++) tmp.push_back(pointsC[i]);
        auto curD = deCast(tmp, t) * double(n);
        tmp.clear();
        rep(i, n) tmp.push_back(pointsC[i]);
        curD -= deCast(tmp, t) * double(n);

        if (norm(res) < eps) {
            if (tt < eps || t < 0 || t > 1) return nullIsect();
            //cout << theta << " " << t <<  " " << a << endl;
            ret.a = 1 / tt;
            ret.hit = r.o + r.dir * tt;
            /* if (Vec3d(ret.hit[0], ret.hit[1], 0).dot(r.dir) > 0) {
                theta = rand() % 1000 * 0.01;
                t = rand() % 1000 * 0.001;
                tt = rand() % 10 * 0.2;
                max_iter = 20;
                continue;
            }*/
            ret.nGeometry = Vec3d(curD[2] * cos(theta), -curD[2] * sin(theta), -curD[0]);
            ret.nGeometry /= norm(ret.nGeometry);
            ret.nShading = ret.nGeometry;
            ret.wo = -r.dir / norm(r.dir);

            ret.isnull = false;
            return ret;
        }


        J(0, 0) = -sin(theta) * cur[0];
        J(1, 0) = cos(theta) * cur[0];
        J(2, 0) = 0;
        J(0, 1) = cos(theta) * curD[0];
        J(1, 1) = sin(theta) * curD[0];
        J(2, 1) = curD[2];
        rep(i, 3) J(i, 2) = -r.dir[i];

        auto dx = -J.inv() * res;

        //cout << J << endl;

        theta += dx[0];
        t += dx[1];
        tt += dx[2];
    }
    return nullIsect();

}