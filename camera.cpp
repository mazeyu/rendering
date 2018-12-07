//
// Created by mzy on 2018/11/5.
//

#include "camera.h"

void pinhole::genRay() {

    int index = 0;
    rep(i, h) {
        rep(j, w) {
            double x = 1 - (double) i / h, y = (double) j / w;
            //x += (rand() % 11 - 5.0) / (40 * h);
            //y += (rand() % 11 - 5.0) / (40 * w);
            sampleX[index] = 1 - x;
            sampleY[index] = y;
            auto dir = (1 - x - y) * axes[0] + x * axes[1] + y * axes[2];
            ray r(o, dir);
            /*rep(ii, 3) {
                rays[index].o.at<double>(ii, 0) = o.at<double>(ii, 0);
                rays[index].dir.at<double>(ii, 0) = dir.at<double>(ii, 0);
            }*/
            rays[index] = r;
            index++;

        }
        cout << i << endl;
    }

}

void camera::synthesize(string filename = "pic") {
    auto n = nRays;
    rep(i, h) rep(j, w) {
        s[i * w + j] = spectrum(0);
        weight[i * w + j] = 0;
    }
    cout << specs[h / 5 * w + w / 2].c[0] << endl;
    double filterR = 0.5;
    rep(i, n) {
        auto x0 = sampleX[i] * h, y0 = sampleY[i] * w;
        auto minx = cvCeil(x0 - filterR), maxx = cvFloor(x0 + filterR),
        miny = cvCeil(y0 - filterR), maxy = cvFloor(y0 + filterR);
        for (int x = minx; x <= maxx; x++)
            for (int y = miny; y <= maxy; y++)
                if (x < h && x >= 0 && y < w && y >= 0) {
                    auto wei = exp(- (x - x0) * (x - x0) - (y - y0) * (y - y0)) - exp(-filterR * filterR);
                    weight[x * w + y] += wei;
                    s[x * w + y] += specs[i] * wei;
                }
    }
    rep(i, h) rep(j, w) s[i * w + j] /= weight[i * w + j];


    rep(i, h) rep(j, w) pic.at<double>(i, j) = s[i * w + j].c[0] * 256;
    imwrite(filename + ".jpg", pic);
}