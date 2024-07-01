# 5-Point Essential Matrix w/ Homotopy Continuation

<!---
TO RUN: Go to 

https://www.unimelb-macaulay2.cloud.edu.au/#tutorial-welcome-1

Upload using "Load your Own Tutorial", and follow instructions to run
code below.
--->

## Goals of this tutorial

1. Illustrate a generic workflow for solving minimal geometric CV problems with Homotopy Continuation (HC).
2. Specifically, we show how to solve the five-point relative pose problem for calibrated cameras.
3. Demonstrate some abilities of Macaulay2, an open-source computer algebra system that can be useful for prototying solvers.

## HC offline / online paradigm

Consider a system of polynomial equations
\[
f_1 (x; p) = \ldots  = f_N (x;p) = 0,
\quad
\underbrace{p \in \mathbb{R}^m}_{\textrm{given}},
\quad 
\underbrace{x \in \mathbb{R}^n}_{\textrm{unknown}}.
\]

We need two technical assumptions:

1. the variety of problem-solution satisfying these equations is irreducible, and
2. the dimension of the variety equals the number of "parameters" m.

In practice, these roughly translate to:
1. it is easy to sample an initial problem-solution pair, and
2. the number of unknown DOF equals the number of input DOF.

Given these two assumptions, HC solvers involve two steps:

1. Offline stage (monodromy): compute all solutions for a generic problem instance $p_0 \in \mathbb{C}^m.$ (This only needs to be done once!)
2. Online stage (parameter homotopy): when presented with a new problem instance $p_1 \in \mathbb{R}^m,$ deform initial solutions of $(f;p_0)=0$ to those of $f(x;p_1)=0.$

Both stages are both implemented in Macaulay2.
Code generation tools (MiNuS by Fabbri, GPU-HC by Chen) allow us to build faster C++ solvers for the online stage.
For the online stage, Macaulay2 can solve systems with 10,000 complex solutions within a few hours.

Selecting "the best" solution is outside the scope of this framework, and thus somebody else's job: eg. scoring in RANSAC, with some cost function, filtering out non-real solutions, neural network, etc.

## Getting started

We load a required Macaulay2 package, which in turn loads several dependancies.

```
restart;
needsPackage "MonodromySolver";
```

Be sure to hit the pulsating START button on the right side of the screen, if you haven't already done so.

## Setting up equations

We start by setting up equations for the five-point problem.
We use the cubic trace constraints, due to Demazure,
\[
2 \mathbf{E} \mathbf{E}^T \mathbf{E} - \operatorname{tr} \left(\mathbf{E} \mathbf{E}^T\right) \mathbf{E} = \mathbf{0}_{3\times 3},
\]
as well as the epipolar constraints for each point correspondence $\mathbf{p_{i 1}} \leftrightarrow \mathbf{p_{i 2}}$,
\[
\mathbf{p_{i 1}}^T \mathbf{E}  \mathbf{p_{i 2}},
\]
plus an additional normalization parameter $s$ to fix the scale of $\mathbf{E}$
\[
s - \left(\mathbf{E_{1,1}} + \cdots + \mathbf{E_{3,3}} \right) = 0.
\]

You can run the code that creates these equations (represented as straight-line programs) here:

```
unknowns = matrix{toList vars(E_(1,1)..E_(3,3))};
params = matrix{toList vars(p_(1,1,1)..p_(5,2,3),s)};
E = gateMatrix for i from 1 to 3 list for j from 1 to 3 list E_(i,j);
EEt = E * transpose E;
demEqs = EEt * E - (1/2) * sum(3, i-> EEt_(i,i)) * E;
ptEqs = for i from 1 to 5 list (matrix{{p_(i,2,1),p_(i,2,2),p_(i,2,3)}} * E * matrix{{p_(i,1,1)},{p_(i,1,2)},{p_(i,1,3)}})_(0,0);
chartEq = s - sum flatten entries E;
eqs = transpose matrix{flatten entries demEqs | ptEqs | {chartEq}};
G = gateSystem(params, unknowns, eqs);
```

## Step 1. (offline) Synthesize an initial problem-solution pair and run monodromy

We now synthesize a random instance of the 5-point problem.
Here is a function that does exactly that.


```
synthesize = FF -> ( -- FF can be real numbers (RR) or complex numbers (CC)
    q0 := apply(5, i -> random(FF^3, FF^1));
    t0 := random(FF^3, FF^1);
    -- can sample rotation w/ Cayley
    (I := id_(FF^3); R0 := random(FF^3, FF^3); R0 = R0 - transpose R0; R0 = (I + R0)*inverse(I-R0));
    p0 := {q0, apply(q0, q -> R0*q+t0)};
    E0 := matrix{{0,-t0_(2,0),t0_(1,0)},{t0_(2,0), 0, -t0_(0,0)},{-t0_(1,0),t0_(0,0),0}} * R0; --; E0 = (1/sum(flatten entries E0)) * E0)
    -- we will normalize E so that the sum of the entries is some fixed constant
    s0 := sum flatten entries E0;
    solution := point{flatten entries E0};
    problem := point{toList apply((1,1,1)..(5,2,3), (i,j,k) -> p0#(j-1)#(i-1)_(k-1,0)) | {s0}};
    (problem, solution)
    )
```

The argument of this function specifies whether we use the real or complex numbers.
For the offline monodromy computation, complex numbers should be used.
Good sanity checks are 

```
(prob0, sol0) = synthesize CC
```

## Step 1 (cont.)

Some good sanity-checking questions before doing any homotopy continuation are:

(1) Does the synthesized solution satisfy all equations?
(2) Are the rank conditions satisfied?

Fortunately, the answer to both questions in this case is yes.

```
(prob0, sol0) = synthesize CC
norm evaluate(G, prob0, sol0)
numericalRank evaluateJacobian(G, prob0, sol0)
```

Once we are satisfied with our ps-pair, it is time to run monodromy.
We select a square subsystem in what follows.
Note, however, that we don't need to worry about extraneous solutions of this subsystem: they are disconnected from the 10 solutions to the true system under the monodromy action.

```
F = squareUp(prob0, sol0, G);
V = first monodromySolve(F, prob0, {sol0});
(prob0, sol0s) = (V.BasePoint, points V.PartialSols)
```

### Step 2. (online) Solve a new problem.

Now let's create a synthetic ps-pair over the real numbers...

```
(prob1, sol1) = synthesize RR
```

and solve it!

```
H01 = specialize(parametricSegmentHomotopy F, transpose(matrix prob0 | matrix prob1))
sol1s = trackHomotopy(H01, sol0s)        
```
Did we recover the ground truth (up to some numerical error)?

```
any(sol1s, x -> norm(matrix x - matrix sol1) < 1e-6)
```
