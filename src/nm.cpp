
/**********************************************************************
    Copyright (c) 2023  Marcel A. Samek. All rights reserved.
    Licensed under the MIT License
    See LICENSE file in the project root for full license information.

    This code was originally derived from the Nelder-Mead algorithm as 
    implemented by Michael F. Hutt whose code can be found at
    https://github.com/huttmf/nelder-mead and was licensed under the
    MIT license and Copyright (c) 1997 Michael F. Hutt
 **********************************************************************/

#include "nm.h"


NelderMead::NelderMead(
    uint32_t inSize,
    const std::function<double(const std::vector<double>&)> & inEvalFunc,
    const std::function<void(std::vector<double>&)> & inConstrainFunc
)
{
    size = inSize;
    evalFunc = inEvalFunc;
    constrainFunc = inConstrainFunc;

    v.resize(size + 1);
    for (uint32_t i = 0; i <= size; i++) {
        v[i].resize(size);
    }
    f.resize(size + 1);
    vr.resize(size);
    ve.resize(size);
    vc.resize(size);
    vm.resize(size);
}

NelderMead::~NelderMead()
{
}


void NelderMead::doInitialize(const std::vector<double> & start, double scale)
{
    double pn, qn;

    evalCount = 0;

    vs = 0;         // vertex with smallest value
    vh = 0;         // vertex with next smallest value
    vg = 0;         // vertex with largest value

    pn = scale * (sqrt(size + 1) - 1 + size) / (size * sqrt(2));
    qn = scale * (sqrt(size + 1) - 1) / (size * sqrt(2));

    for (uint32_t i = 0; i < size; i++) {
        v[0][i] = start[i];
    }

    for (uint32_t i = 1; i <= size; i++) {
        for (uint32_t j = 0; j < size; j++) {
            if (i - 1 == j) {
                v[i][j] = pn + start[j];
            }
            else {
                v[i][j] = qn + start[j];
            }
        }
    }
}


#if NELDER_MEAD_DEBUG
void NelderMead::doPrintStart()
{
    printf("Initial Values\n");
    for (uint32_t j = 0; j <= size; j++) {
        for (uint32_t i = 0; i < size; i++) {
            printf("%f, ", v[j][i]);
        }
        printf("value %f\n", f[j]);
    }
}

void NelderMead::doPrintIteration(uint32_t itr)
{
    printf("Iteration %d\n", itr);
    for (uint32_t j = 0; j <= size; j++) {
        for (uint32_t i = 0; i < size; i++) {
            printf("%f %f\n", v[j][i], f[j]);
        }
    }
}
#endif

void NelderMead::doIndexes()
{
    vh = vs;
    for (uint32_t j = 0; j <= size; j++) {
        if (f[j] > f[vg]) {
            vg = j;
        }
        if (f[j] < f[vs]) {
            vs = j;
        }
    }
    vh = vs;
    for (uint32_t j = 0; j <= size; j++) {
        if (f[j] > f[vh] && f[j] < f[vg]) {
            vh = j;
        }
    }
}

void NelderMead::exec(const std::vector<double> & start, double tolerancee, double scale)
{
    double fr;      // value of function at reflection point
    double fe;      // value of function at expansion point
    double fc;      // value of function at contraction point

    // This function can be called many times for the same instance of the class
    // so we have to initialize it every time.
    doInitialize(start, scale);

    // The starting values that we were passed might not actually obey the constraint
    // function that was provided. Silly caller. So we constrain them here.
    if (constrainFunc) {
        for (uint32_t j = 0; j <= size; j++) {
            constrainFunc(v[j]);
        }
    }

    // find the initial function values based on the freshly constraine starting values
    for (uint32_t j = 0; j <= size; j++) {
        f[j] = doEvaluate(v[j]);
    }

#if NELDER_MEAD_DEBUG
    // print out the initial values
    doPrintStart();
#endif

    // The loop that converges (maybe) on a what is being sought
    uint32_t iterationCount = 0;
    while (++iterationCount <= configMaxIterations) {

        // calculate the  indexes of significant vertices of the simplex
        // that will be used in subsequent calculations. 
        doIndexes();

        // calculate the centroid of the simplex
        double cent;
        for (uint32_t j = 0; j <= size - 1; j++) {
            cent = 0.0;
            for (uint32_t m = 0; m <= size; m++) {
                if (m != vg) {
                    cent += v[m][j];
                }
            }
            vm[j] = cent / size;
        }

        // reflect vg to new vertex vr. The reflection might need to be constrained.
        for (uint32_t j = 0; j <= size - 1; j++) {
            vr[j] = vm[j] + configReflectionCoefficient * (vm[j] - v[vg][j]);
        }
        if (constrainFunc) {
            constrainFunc(vr);
        }

        // recalculate the simplex values
        fr = doEvaluate(vr);
        if (fr < f[vh] && fr >= f[vs]) {
            for (uint32_t j = 0; j <= size - 1; j++) {
                v[vg][j] = vr[j];
            }
            f[vg] = fr;
        }

        // investigate a step further in this direction 
        if (fr < f[vs]) {
            for (uint32_t j = 0; j <= size - 1; j++) {
                ve[j] = vm[j] + configExpansionCoefficient * (vr[j] - vm[j]);
            }
            if (constrainFunc != NULL) {
                constrainFunc(ve);
            }
            fe = doEvaluate(ve);

            if (fe < fr) {
                for (uint32_t j = 0; j <= size - 1; j++) {
                    v[vg][j] = ve[j];
                }
                f[vg] = fe;
            }
            else {
                for (uint32_t j = 0; j <= size - 1; j++) {
                    v[vg][j] = vr[j];
                }
                f[vg] = fr;
            }
        }

        // check to see if a contraction is necessary 
        if (fr >= f[vh]) {
            if (fr < f[vg] && fr >= f[vh]) {
                // perform outside contraction 
                for (uint32_t j = 0; j <= size - 1; j++) {
                    vc[j] = vm[j] + configContractionCoefficient * (vr[j] - vm[j]);
                }
                if (constrainFunc != NULL) {
                    constrainFunc(vc);
                }
                fc = doEvaluate(vc);
            }
            else {
                // perform inside contraction 
                for (uint32_t j = 0; j <= size - 1; j++) {
                    vc[j] = vm[j] - configContractionCoefficient * (vm[j] - v[vg][j]);
                }
                if (constrainFunc != NULL) {
                    constrainFunc(vc);
                }
                fc = doEvaluate(vc);
            }


            if (fc < f[vg]) {
                for (uint32_t j = 0; j <= size - 1; j++) {
                    v[vg][j] = vc[j];
                }
                f[vg] = fc;
            }

            else {
                // at this point the contraction is not successful,
                // we must halve the distance from vs to all the
                // vertices of the simplex and then continue.
                for (uint32_t row = 0; row <= size; row++) {
                    if (row != vs) {
                        for (uint32_t j = 0; j <= size - 1; j++) {
                            v[row][j] = v[vs][j] + (v[row][j] - v[vs][j]) / 2.0;
                        }
                    }
                }

                // re-evaluate all the vertices 
                for (uint32_t j = 0; j <= size; j++) {
                    f[j] = doEvaluate(v[j]);
                }

                // calculate significant indexes of the simplex
                doIndexes();

                if (constrainFunc != NULL) {
                    constrainFunc(v[vg]);
                }
                f[vg] = doEvaluate(v[vg]);
                if (constrainFunc != NULL) {
                    constrainFunc(v[vh]);
                }
                f[vh] = doEvaluate(v[vh]);
            }
        }

        // print out the value at each iteration
#if NELDER_MEAD_DEBUG
        doPrintIteration(iterationCount);
#endif

        // test for convergence
        double fsum = 0.0;
        for (uint32_t j = 0; j <= size; j++) {
            fsum += f[j];
        }
        double favg = fsum / (size + 1);
        double s = 0.0;
        for (uint32_t j = 0; j <= size; j++) {
            s += pow((f[j] - favg), 2.0) / (size);
        }
        s = sqrt(s);

        if (s < tolerancee) {
            break;
        }
    }

    // calculate significant indexes of the simplex
    doIndexes();

    // evaluate the minimum  and stuff the results
    lastExecResults.min = doEvaluate(v[vs]);
    lastExecResults.evalCount = evalCount;
    lastExecResults.iterationCount = iterationCount;
    lastExecResults.minValues.clear();
    for (uint32_t j = 0; j < size; j++) {
        lastExecResults.minValues.push_back(v[vs][j]);
    }
}

