
/**********************************************************************
    Copyright (c) 2023  Marcel A. Samek. All rights reserved.
    Licensed under the MIT License
    See LICENSE file in the project root for full license information.
 **********************************************************************/

#include "nm.h"

#include <math.h>
#include <stdio.h>


double myFunction(const std::vector<double> & x)
{
    double a = -1.23456;
    double b = 6.54321;

    double v1 = b * b - a;
    double v2 = x[0] * x[0] - x[1];

    double w1 = a * a - b;
    double w2 = x[1] * x[1] - x[0];

    double v = sqrt((v1-v2)*(v1-v2) + (w1 - w2) * (w1 - w2));

    return v;
}
void myConstraints(std::vector<double>& x)
{
    for (auto & xi : x) {
        if (xi < -600) {
            xi = -600;
        }
        if (xi > 600) {
            xi = 600;
        }
    }
}

void printResults(const NelderMeadResults& results)
{
    printf("    %u Function Evaluations\n", results.evalCount);
    printf("    %u Iterations through program\n", results.iterationCount);
    printf("    Best result: %le\n", results.min);
    for (uint32_t i = 0; i < results.minValues.size(); i++) {
        printf("        Best variables: %le\n", results.minValues[i]);
    }
}

int main()
{
    // allocate the Nelder-Mead object searching for 2 variables
    NelderMead* simp = new NelderMead(2, myFunction, nullptr);

    simp->setMaxIterations(100000);

    fprintf(stdout, "Trying Nelder Mead with tolerance 1.0e-6\n");
    simp->exec(std::vector<double>{ 1, 1 }, 1.0e-6, 1.0);
    printResults(simp->getLastExecResults());

    fprintf(stdout, "Trying Nelder Mead with tolerance 1.0e-12\n");
    simp->exec(std::vector<double>{ 1, 1 }, 1.0e-12, 1.0);
    printResults(simp->getLastExecResults());

    return 0;
}
