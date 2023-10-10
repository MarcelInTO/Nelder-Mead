
/**********************************************************************
    Copyright (c) 2023  Marcel A. Samek. All rights reserved.
    Licensed under the MIT License
    See LICENSE file in the project root for full license information.

    This code was originally derived from the Nelder-Mead algorithm as
    implemented by Michael F. Hutt whose code can be found at
    https://github.com/huttmf/nelder-mead and was licensed under the
    MIT license and Copyright (c) 1997 Michael F. Hutt
 **********************************************************************/

#pragma once

// system headers
#include <stdint.h>

// std library headers
#include <functional>
#include <vector>

// Set to 1 to enable debug output
#define NELDER_MEAD_DEBUG 0


// Holds the results of the last completed exec call
struct NelderMeadResults {
    uint32_t iterationCount;
    uint32_t evalCount;
    std::vector<double> minValues;
    double min;
};

// Main interface for the Nelder-Mead algorithm. Once constructed, the number of variables
// being solved for cannot be changed.
class NelderMead
{
    public:
        // Constructors and destructor

        NelderMead(
            uint32_t inSize,
            const std::function<double(const std::vector<double>&)> &,
            const std::function<void(std::vector<double>&)> &
        );
        ~NelderMead();

        // public methods

        void exec(const std::vector<double> & inStart, double tolerance, double scale);
        const NelderMeadResults & getLastExecResults() const { return lastExecResults; }
        void setMaxIterations(uint32_t inValue) { configMaxIterations = inValue; }
        void setReflectionCoefficient(double inValue) { configReflectionCoefficient = inValue; }
        void setContractionCoefficient(double inValue) { configContractionCoefficient = inValue; }
        void setExpansionCoefficient(double inValue) { configExpansionCoefficient = inValue; }


    private:
        // Configuration values that can be modified by the user prior to an exec call

        uint32_t configMaxIterations = 1000;
        double configReflectionCoefficient = 1.0;
        double configContractionCoefficient = 0.5;
        double configExpansionCoefficient = 2.0;

        // Core definition of an instantiation of the algorithm. These
        // values are set at construction time and cannot be modified 
        // once set. If they need to change, a new instance of the class
        // should be allocated with appropriate values

        uint32_t size = 0;
        std::function<double(const std::vector<double>&)> evalFunc;
        std::function<void(std::vector<double>&)> constrainFunc;

        std::vector<std::vector<double>> v;     // holds vertices of simplex 
        std::vector<double> f;      // value of function at each vertex 
        std::vector<double> vr;     // reflection - coordinates 
        std::vector<double> ve;     // expansion - coordinates 
        std::vector<double> vc;     // contraction - coordinates 
        std::vector<double> vm;     // centroid - coordinates 

        // Current execution state. Reset on every exec call.

        mutable uint32_t evalCount = 0;
        uint32_t vs = 0;         // index of vertex with smallest value
        uint32_t vh = 0;         // index of vertex with next smallest value
        uint32_t vg = 0;         // index of vertex with largest value 

        NelderMeadResults lastExecResults;

        // private methods

        double doEvaluate(const std::vector<double>& x) const
        {
            evalCount++;
            return evalFunc(x);
        }

        void doInitialize(const std::vector<double>& start, double scale);
        void doCentroid(uint32_t vg);
        void doIndexes();

    #if NELDER_MEAD_DEBUG
        void doPrintStart();
        void doPrintIteration(uint32_t itr);
    #endif
};
