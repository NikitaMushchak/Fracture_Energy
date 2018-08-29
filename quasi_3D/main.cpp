#include <stdlib.h>
#include "util.h"

//void Icos1();
//void Modeling3D();
////double Modeling3D_1(double);
//void Modeling3D_2();
//void Modeling3D_3();
//void Modeling3D_1_1();
//void Modeling3D_1_2();
//void Modeling3D_4();
//void Modeling3D_5();
//void Modeling3D_6();
//void Modeling3D_7();
void Modeling3Dv2();
//void ell_test();
//void test();

int main()
{
    if (system("cd ./output"))
    {
        system("mkdir output");
        print("Создаем папку \"output\"\n");

        if (system("cd ./output"))
            print("Папка \"output\" не видна\n");
    }

    //Icos1();
    //Modeling3D();
    //Modeling3D_1(0.9);
    //Modeling3D_1_2();
    //Modeling3D_7();
    Modeling3Dv2();
    //ell_test();
    //test();
    //system("pause");
    return 0;
}
