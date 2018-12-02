#include <iostream>
#include "core.h"
#include "scene.h"





int main() {
    /*vector<Vec3d> v{Vec3d(0, 0, 0), Vec3d(1, 0, 0), Vec3d(1, 0, 1),
                    Vec3d(0, 0, 1)};
    ray r(Vec3d(1.5, 0, -1.0001), Vec3d(-1, 0, 0));

    revoBezier test(v);
    double a;
    cout << test.intersect(r)->a << endl;*/

    /*scene scene1;
    scene1.init();
    scene1.buildKdTree();
    scene1.SPPM();*/
    int sum = 0;
    for (int i = 0; i < 100000000; i++)
        for (int j = 0; j < 10; j++)
            sum += i * j;
}