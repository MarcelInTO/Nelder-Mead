# Nelder-Mead C++
A basic implementation of the Nelder Mead algorightm in C++.

## Overview

Using this implementaiton of Nelder Mead is a two step process:

1. Allocate a NelderMeade solver object
2. Execute the minimization search on that object.

The first step allocates a solver instance that has all the necessary memory allocated.

By separating the actual minimization search into a separate step, we allow the caller to quickly redo the search using different starting points and execution parameters without having to perform any memory allocations or reallocations.

In order to allocate and use a solver, the caller must provide a pointer to an evaluation function. This pointer is provied as a std::function object so that it can either be a plain function or a class method with an instance pointer.

The call to execute the solver allows the caller to specify the starting point for the search, a scale to use for setting the size of the initial simplex, as well as the tolerance for the search

```
        void exec(
            const std::vector<double> & startValues, 
            double tolerance, 
            double scale
        );
```

The solver object provides methods to change several execution parameters after it has been allocated. This allows the caller to call exec() many times, using different settings, without having to reallocate memory.

* Maximum number of iterations
* Reflection coefficient
* Contraction coefficient
* Expansion Coefficient

## Minimal Example

Here is a simple example which allocates a solver and then calls it twice, each time with a different tolerance, so we can see the difference between the number of iterations it took for each tolerance value.

Included in the project is a slightly more complete example.

```
#include "nm.h"

#include <stdio.h>

int main()
{
    // allocate the Nelder-Mead object searching for 2 variables
    NelderMead* simp = new NelderMead(2, myFunction, nullptr);

    // set the number of iterations to a non-default number
    simp->setMaxIterations(500);

    // execute the search a print the result
    simp->exec(std::vector<double>{ 1, 1 }, 1.0e-6, 1.0);
    printResults(simp->getLastExecResults());

    // sececute the search with a different tolerance
    simp->exec(std::vector<double>{ 1, 1 }, 1.0e-12, 1.0);
    printResults(simp->getLastExecResults());

    return 0;

}
```

Note that the results of an exec() call on the solver are retrieved by calling the getLastExectResult() method. In this example there is a simple function that prints the output of those results:


```
void printResults(const NelderMeadResults& results)
{
    printf("    %u Function Evaluations\n", results.evalCount);
    printf("    %u Iterations through program\n", results.iterationCount);
    printf("    Best result: %le\n", results.min);
    for (uint32_t i = 0; i < results.minValues.size(); i++) {
        printf("        Best variables: %le\n", results.minValues[i]);
    }
}
```

In order to allocate the solver we had to pass it our evaluation function. 

```
// my test function with two variables
double myFunction(const std::vector<double> & x)
{
    ..My code to evaluate function...

    return v;
}
```


