Selected code from: [Murphy, Robertson, Harty & O'Connell (2015). Neural evidence accumulation persists after choice to inform metacognitive judgments, _eLife_, 4: e11946](https://elifesciences.org/articles/11946#abstract).

### Contents

* 1choice_DDM: Matlab functions to fit the one-choice drift diffusion model to error detection/response time data (e.g. as in [Figure 6](https://elifesciences.org/articles/11946/figures) in paper). Top-level function is `FIT_diffusion_model_PSO_SS.m`. Optimization performed via particle swarm optimization ([Birge, 2003](https://ieeexplore.ieee.org/document/1202265)). Model predictions for a given parameter set generated via Monte Carlo simulation, since one-choice DDM has no known analytical solution. Code includes optional accumulator leak, bound collapse and pre-error accumulation onset parameters, which I played around with at the time but did not include in the paper for various reasons. N.B. simulation code is not at all optimized - this was my first effort at modelling and it's very inefficient :-D A small amount of effort will speed up simulations significantly.

### License
Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions: The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software. If you use the Software for your own research, cite the paper.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
