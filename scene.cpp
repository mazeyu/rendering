//
// Created by mzy on 2018/11/5.
//

#include "scene.h"
#include "BSDF.h"

double randomReal() {
    return rand() % 100000 * 0.00001;
}

Vec3d sphereSampler() {
    double x, y, z;
    for (;;) {
        x = randomReal() * 2 - 1;
        y = randomReal() * 2 - 1;
        z = randomReal() * 2 - 1;
        if (x * x + y * y + z * z < 1) {
            return Vec3d(x, y, z) / norm(Vec3d(x, y, z));
        }
    }

}

struct SPPMPixel {
    struct VisiblePoint {
        VisiblePoint() {}

        VisiblePoint(Vec3d p, Vec3d wo, BSDF *bsdf, spectrum beta) :
                p(p), wo(wo), bsdf(bsdf), beta(beta) {}

        Vec3d p, wo;
        const BSDF *bsdf = nullptr;
        spectrum beta;
    } vp;

    SPPMPixel() : M(0) {}

    atomic<int> M;
    spectrum tau, Ld, phi;
    double N = 0, radius = 0;
};

struct SPPMPixelListNode {
    SPPMPixel *pixel;
    SPPMPixelListNode *next;
};


void free(SPPMPixelListNode *&cur) {
    if (cur == nullptr) return;
    free(cur->next);
    delete cur;
    cur = nullptr;
}

inline unsigned int hashP(const Vec3i &p, int hashSize) {
    return (unsigned int) ((p[0] * 73856093) ^ (p[1] * 19349663) ^
                           (p[2] * 83492791)) %
           hashSize;
}

spectrum scene::UniformSampleOneLight(intersection &isect) {
    auto n = lights.size(), k = rand() % n;
    Vec3d wi;
    double pdf;
    auto Li = lights[k]->sampleLi(isect, wi, pdf);
    auto isect1 = intersect(ray(isect.hit, wi));
    if (!isect1.isnull && isect1.a > 1) return spectrum();
    wi /= norm(wi);
    auto L = Li * isect.bsdf->f(isect.wo, wi) * abs(wi.dot(isect.nShading)) / pdf * n;
    return L;

}

void scene::SPPM() {
    cout << "SPPMing" << endl;
    freopen("log.txt", "w", stdout);
    int iteration = 50, photonsPerIter = 10000, maxDep = 10;

    cam->genRay();

    auto rays = cam->rays;
    auto n = cam->nRays;
    SPPMPixel *pixels = new SPPMPixel[n];

    rep(i, n) pixels[i].radius = 0.1;//(kdRoot->BB.max[0] - kdRoot->BB.min[0]) / 20;

    int nThread = 1;
    thread t[nThread];
    auto batch = (n + nThread - 1) / nThread;
    rep(iter, iteration) {
        spectrum beta(1);
        rep(i, nThread) {
            t[i] = thread([&, i]() {
                rep(j, batch)
                    if (i * batch + j < n) {
                        if (j % 100 == 0) cout << j << endl;
                        int index = i * batch + j;
                        auto r = rays[index];
                        cout << "l1" << endl;
                        for (int depth = 0; depth < maxDep; depth++) {
                            cout << "l2" << endl;
                            auto isect = intersect(r);
                            cout << "l3" << endl;
                            if (isect.isnull) {
                                for (const auto &light : lights)
                                    pixels[index].Ld += light->Le(r) * beta;
                                break;
                            }
                            cout << "l4" << endl;
                            const BSDF &bsdf = *isect.bsdf;
                            auto wo = -r.dir;
                            pixels[index].Ld += beta * UniformSampleOneLight(isect);

                            cout << "l5" << endl;
                            auto isDiffuse = true;

                            cout << "l6" << endl;
                            if (isDiffuse) {
                                pixels[index].vp.p = isect.hit;
                                pixels[index].vp.wo = wo;
                                pixels[index].vp.bsdf = &bsdf;
                                pixels[index].vp.beta = beta;
                                cout << "l7" << endl;
                                break;
                            }
                            ////////////////////

                        }

                    }
            });
        }
        rep(i, nThread) t[i].join();


        int hashSize = n;
        vector<atomic<SPPMPixelListNode *> > grid(hashSize);

        //int gridRes[3];
        double maxRadius = 0;
        rep(i, n) maxRadius = max(maxRadius, pixels[i].radius);
        //rep(i, 3) gridRes[i] = cvCeil((kdRoot->BB.max[i] - kdRoot->BB.min[i]) / maxRadius);

        rep(i, nThread) {
            t[i] = thread([&, i]() {
                rep(j, batch)if (i * batch + j < n) {
                        int index = i * batch + j;
                        SPPMPixel &pixel = pixels[index];
                        auto radius = pixel.radius;
                        Vec3i pMin, pMax;
                        // because of hash, no overflow
                        rep(k, 3) {
                            pMin[k] = cvFloor((pixel.vp.p[k] - radius) / maxRadius);
                            pMax[k] = cvFloor((pixel.vp.p[k] + radius) / maxRadius);
                        }
                        replr(x, pMin[0], pMax[0])
                            replr(y, pMin[1], pMax[1])
                                replr(z, pMin[2], pMax[2]) {
                                    int h = hashP(Vec3i(x, y, z), hashSize);
                                    auto *node = new SPPMPixelListNode;
                                    node->pixel = &pixel;
                                    while (!grid[h].compare_exchange_weak(
                                            node->next, node));
                                }

                    }
            });
        }
        rep(i, nThread) t[i].join();

        auto batchP = (photonsPerIter + nThread - 1) / nThread;

        cout << endl << "emitting photon" << endl;
        rep(i, nThread) {
            t[i] = thread([&, i]() {
                rep(j, batchP)if (i * batchP + j < photonsPerIter) {

                        if (j % 50 == 0) cout << j << " ";
                        auto nL = lights.size(), k = rand() % nL;
                        ray r;
                        Vec3d normal;
                        double pdfDir, pdfPos;
                        auto Le = lights[k]->sampleLe(r, normal, pdfDir, pdfPos);
                        auto beta = Le * abs(normal.dot(r.dir)) * nL / (pdfPos * pdfDir);

                        for (int depth = 0; depth < maxDep; depth++) {
                            auto isect = intersect(r);
                            if (isect.isnull) break;
                            if (depth > 0) {
                                Vec3i posG;
                                rep(d, 3) posG[d] = cvFloor(isect.hit[d] / maxRadius);
                                auto h = hashP(posG, hashSize);
                                for (SPPMPixelListNode *node = grid[h].load(std::memory_order_relaxed);
                                     node != nullptr; node = node->next) {
                                    SPPMPixel &pixel = *node->pixel;
                                    auto radius = pixel.radius;
                                    if (norm(pixel.vp.p - isect.hit) > radius) continue;
                                    Vec3d wi = -r.dir;
                                    auto phi = beta * pixel.vp.bsdf->f(pixel.vp.wo, wi);
                                    pixel.phi += phi;
                                    ++pixel.M;
                                }
                            }
                            Vec3d wi, wo = -r.dir;
                            double pdf;
                            auto fr = isect.bsdf->sampleF(wo, wi, pdf, isect.nShading);
                            auto bnew = beta * fr * abs(wi.dot(isect.nShading)) / pdf;
                            double q = max(0., 1 - bnew.y() / beta.y());
                            if (randomReal() < q) break;
                            beta = bnew / (1 - q);
                            r = ray(isect.hit, wi);


                        }

                    }
            });
        }
        rep(i, nThread) t[i].join();
        rep(i, hashSize) free(grid[i]);
        rep(i, n) {
            SPPMPixel &p = pixels[i];
            if (p.M > 0) {
                double gamma = 2.0 / 3, Nnew = p.N + gamma * p.M;
                auto Rnew = p.radius * sqrt(Nnew / (p.N + p.M));
                auto phi = p.phi;
                p.tau = (p.tau + p.vp.beta * phi) * (Rnew * Rnew) /
                        (p.radius * p.radius);

                p.N = Nnew;
                p.radius = Rnew;
                p.M = 0;
                p.phi = spectrum();
            }
            p.vp.beta = spectrum();
            p.vp.bsdf = nullptr;
        }


        int writeFreq = 1;
        int64 Np = (iter + 1) * photonsPerIter;
        if (iter % writeFreq == 0) {
            rep(i, n) {
                auto L = pixels[i].Ld;
                L /= iter + 1;
                L += pixels[i].tau / (Np * acos(-1) * pixels[i].radius * pixels[i].radius);

                cam->specs[i] = L;
            }
            cam->synthesize(to_string(iter));
        }

    }
}


void scene::render() {
    cout << 0 << endl;

    cam->genRay();
    auto rays = cam->rays;
    auto n = cam->nRays;
    cout << 1 << endl;

    int nThread = 100;
    thread t[nThread];
    auto batch = (n + nThread - 1) / nThread;
    rep(i, nThread) {
        t[i] = thread([&, i]() {
            rep(j, batch)if (i * batch + j < n) {
                    if (j % 100 == 0) cout << j << endl;
                    cam->specs[i * batch + j] = integL(rays[i * batch + j]);
                }
        });
    }
    rep(i, nThread) t[i].join();
    cout << 2 << endl;
    cam->synthesize("rt");

}

spectrum scene::RT(ray r) {
    spectrum ret;
    auto intersection = intersect(r);
    if (intersection.isnull) return ret;
    auto hit = intersection.hit;
    auto n = intersection.nShading;
    auto spec = intersection.spec;

    ret = UniformSampleOneLight(intersection);
    return ret;
}

spectrum scene::integL(ray r) {

    spectrum s;
    s = RT(r);
    //rep(ii, 3) cout << r.dir.at<double>(ii, 0) << " ";
    //cout << endl;
    return s;
}


void scene::addBezier() {
    ifstream in("model/mushroom.txt");
    string s;
    vector<Vec3d> points;
    while (getline(in, s)) {
        stringstream ss(s);
        double x, y, z;
        ss >> x >> y >> z;
        z = -z;
        points.emplace_back(Vec3d(x, y, z));
    }
    rep(i, (points.size() - 1) / 3 + 1) {
        vector<Vec3d> tmp{points[i * 3], points[i * 3 + 1], points[i * 3 + 2], points[i * 3 + 3]};
        int divTheta = 10, divT = 2;
        rep(iTheta, divTheta) rep(iT, divT) rep(k, 2)objs.emplace_back(
                            new primitive(subRevoBezier(tmp, iTheta, divTheta, iT, divT, k)));
    }
}

void scene::init() {
    kdRoot = nullptr;
    double o[] = {3, 2, 1};
    double at[] = {0, 0, 0};
    double up[] = {0, 1, 0};
    cam = new pinhole(Vec3d(o), Vec3d(at), Vec3d(up));

    //lights.emplace_back(new pointLight(Vec3d(20, 0, 0)));
    lights.emplace_back(new pointLight(Vec3d(3, 2, 1), 5));

    //addBezier();
    addFromFile();
    /*vector<Vec3d> v{Vec3d(0, 0, 0), Vec3d(1, 0, 0), Vec3d(1, 0, 1),
                    Vec3d(0, 0, 1)};


    objs.emplace_back(new primitive(new revoBezier(v)));*/

    int sz = 5, ht =  5;

    objs.emplace_back(new primitive(new triangle(Vec3d(-sz, -ht, -sz),
                                                 Vec3d(-sz, -ht, sz),
                                                 Vec3d(sz, -ht, -sz))));


    objs.emplace_back(new primitive(new triangle(Vec3d(sz, -ht, sz),
                                                 Vec3d(sz, -ht, -sz),
                                                 Vec3d(-sz, -ht, sz))));

    objs.emplace_back(new primitive(new triangle(Vec3d(-sz, ht, -sz),
                                                 Vec3d(-sz, ht, sz),
                                                 Vec3d(sz, ht, -sz))));


    objs.emplace_back(new primitive(new triangle(Vec3d(sz, ht, sz),
                                                 Vec3d(sz, ht, -sz),
                                                 Vec3d(-sz, ht, sz))));

    objs.emplace_back(new primitive(new triangle(Vec3d(-ht, -sz, -sz),
                                                 Vec3d(-ht, -sz, sz),
                                                 Vec3d(-ht, sz, -sz))));


    objs.emplace_back(new primitive(new triangle(Vec3d(-ht, sz, sz),
                                                 Vec3d(-ht, sz, -sz),
                                                 Vec3d(-ht, -sz, sz))));

    objs.emplace_back(new primitive(new triangle(Vec3d(ht, -sz, -sz),
                                                 Vec3d(ht, -sz, sz),
                                                 Vec3d(ht, sz, -sz))));


    objs.emplace_back(new primitive(new triangle(Vec3d(ht, sz, sz),
                                                 Vec3d(ht, sz, -sz),
                                                 Vec3d(ht, -sz, sz))));

    objs.emplace_back(new primitive(new triangle(Vec3d(-sz, -sz, -ht),
                                                 Vec3d(-sz, sz, -ht),
                                                 Vec3d(sz, -sz, -ht))));


    objs.emplace_back(new primitive(new triangle(Vec3d(sz, sz, -ht),
                                                 Vec3d(sz, -sz, -ht),
                                                 Vec3d(-sz, sz, -ht))));

    objs.emplace_back(new primitive(new triangle(Vec3d(-sz, -sz, ht),
                                                 Vec3d(-sz, sz, ht),
                                                 Vec3d(sz, -sz, ht))));


    objs.emplace_back(new primitive(new triangle(Vec3d(sz, sz, ht),
                                                 Vec3d(sz, -sz, ht),
                                                 Vec3d(-sz, sz, ht))));

}

void scene::addFromFile() {
    ifstream in("model/pigD.obj");
    string s;
    vector<Vec3d> points, nPoints;

    while (getline(in, s)) {
        if (s.substr(0, 2) == "v ") {
            double x[3];
            sscanf(s.c_str(), "v %lf%lf%lf", &x[0], &x[1], &x[2]);
            points.emplace_back(Vec3d(x));
        }
        if (s.substr(0, 3) == "vn ") {
            double x[3];
            sscanf(s.c_str(), "vn %lf%lf%lf", &x[0], &x[1], &x[2]);
            nPoints.emplace_back(Vec3d(x));
        }
        if (s.substr(0, 2) == "f ") {
            string t;
            stringstream ss(s);
            vector<int> inds, inds2;
            while (ss >> t) {
                if (t != "f") {
                    int ind, ind2, ind3;
                    sscanf(t.c_str(), "%d/%d/%d", &ind, &ind3, &ind2);
                    inds.emplace_back(ind - 1);
                    inds2.emplace_back(ind2 - 1);
                }
            }

            rep(i, inds.size() - 2) {
                Vec3d tmp[3]{points[inds[0]], points[inds[i + 1]], points[inds[i + 2]]};
                Vec3d ntmp[3]{nPoints[inds2[0]], nPoints[inds2[i + 1]], nPoints[inds2[i + 2]]};
                objs.emplace_back(new primitive(new triangle(tmp, ntmp)));
            }
        }
    }
    cout << objs.size() << endl;

}

bool bbIntersect(bb BB, ray r) {
    double mind = 0, maxd = oo;
    rep(i, 3) {
        double d1 = (BB.min[i] - r.o[i]) / r.dir[i];
        double d2 = (BB.max[i] - r.o[i]) / r.dir[i];
        if (d1 > d2) swap(d2, d1);
        if (d1 > mind) mind = d1;
        if (d2 < maxd) maxd = d2;
        if (mind > maxd - eps) return false;
    }
    return true;
}

vector<primitive *> scene::kdIntersect(kdNode *cur, ray r) {
    vector<primitive *> ret;
    if (!bbIntersect(cur->BB, r)) return ret;
    if (cur->lc == nullptr) return cur->prims;
    auto pl = kdIntersect(cur->lc, r), pr = kdIntersect(cur->rc, r);
    ret.insert(ret.end(), pl.begin(), pl.end());
    ret.insert(ret.end(), pr.begin(), pr.end());
    return ret;
}

intersection scene::intersect(ray r) {
    vector<primitive *> toInter;
    if (kdRoot == nullptr) toInter = objs;
    else toInter = kdIntersect(kdRoot, r);
    priority_queue<pair<double, intersection> > candidate;
    rep(i, toInter.size()) {
        auto ret = toInter[i]->Shape->intersect(r);
        if (!ret.isnull)
            candidate.push(mp(ret.a, ret));
    }
    if (candidate.empty()) return nullIsect();
    return candidate.top().second;
}

void build(kdNode *&cur, bb BB, vector<primitive *> prims, int &kdCnt) {
    kdCnt += prims.size();
    cur = new kdNode(BB, prims);
    cur->lc = cur->rc = nullptr;
    int maxDim;
    minMaxIdx(BB.max - BB.min, nullptr, nullptr, nullptr, &maxDim);
    //print(BB.max-BB.min);
    //cout << maxDim;
    double pos = (BB.max[maxDim] + BB.min[maxDim]) / 2;
    auto BBl = BB, BBr = BB;
    BBl.max[maxDim] = pos;
    BBr.min[maxDim] = pos;
    vector<primitive *> pl, pr;
    rep(i, prims.size()) {
        auto flag = prims[i]->Shape->intersectPlane(maxDim, pos);
        if (flag & 1) pl.emplace_back(prims[i]);
        if (flag & 2) pr.emplace_back(prims[i]);
    }
    if (pl.size() == prims.size() || prims.size() == pr.size()) return;

    cout << prims.size() << " " << maxDim << pl.size() << ":" << pr.size() << endl;
    //print(BB.min);
    //print(BB.max);

    if (pl.size()) build(cur->lc, BBl, pl, kdCnt);
    if (pr.size()) build(cur->rc, BBr, pr, kdCnt);
}

void scene::buildKdTree() {
    Mat all(3, objs.size() * 2, CV_64FC1);
    double minx[3], maxx[3];
    rep(i, objs.size())
        rep(j, 3) {
            all.at<double>(j, i * 2) = objs[i]->Shape->minmax()(j, 0);
            all.at<double>(j, i * 2 + 1) = objs[i]->Shape->minmax()(j, 1);
        }
    //print(all);
    rep(i, 3) minMaxIdx(all.row(i), &minx[i], &maxx[i]);
    auto BB = bb(Vec3d(minx), Vec3d(maxx));
    print(BB.min);
    print(BB.max);
    build(kdRoot, BB, objs, kdCnt);
    cout << kdCnt << endl;
}