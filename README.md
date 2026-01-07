# Temporal Mapper Tutorial: Lotka-Volterra System

This folder contains a simplified example of temporal mapper construction using the classic **Lotka-Volterra predator-prey model**. The goal is to understand the tmapper methodology before applying it to more complex neural data.

Lotka-Volterra has simple dynamics (i.e. two state variables), a clear periodic structure (i.e. limit cycle behavior) and it is easy to visualize and verify results.


## The Lotka-Volterra Model

The non-linear ODE model is defined as:

$$
\dfrac{\mathrm{d}x}{\mathrm{d}t} = \alpha x - \beta xy \quad \text{(prey: growth minus predation)}

\dfrac{\mathrm{d}y}{\mathrm{d}t} = \delta xy - \gamma y \quad \text{(predator: growth from eating minus death)}
$$


Where:
- $x$ or $x(t)$ denotes prey population
- $y$ or $y(t)$ denotes predator population
- $\alpha$ represents prey birth rate
- $\beta$ represents predation rate
- $\delta$ represents predator birth rate (from consuming prey)
- $\gamma$ represents predator death rate

With default parameters, the system exhibits **limit cycle oscillations** around a fixed point.

## Files in This Folder

| File | Description |
|------|-------------|
| `LotkaVolterra.m` | MATLAB class for simulating the model |
| `run_lotka_volterra.m` | Main tutorial script - **start here!** |
| `parameter_sweep.m` | Explores effect of different parameters |
| `plot_mapper_results.m` | Visualization utility |
| `plot_shape_graph_animated.m` | Animation of trajectory on shape graph |

## How to navigate the repo?

1. Open MATLAB and navigate to this folder
2. Run the main tutorial:
   ```matlab
   run_lotka_volterra
   ```
3. Explore parameter effects:
   ```matlab
   parameter_sweep
   ```

## Understanding Temporal Mapper

### The Key Steps

1. **Simulate/obtain time series data**
   - In our case: prey and predator populations over time

2. **Construct temporal k-NN graph** (`tknndigraph`)
   - Connect each point to its `k` nearest neighbors in state space
   - Also connect consecutive time points (temporal neighbors)
   - This captures both geometry AND temporal ordering

3. **Filter the graph** (`filtergraph`)
   - Points within distance `d` are collapsed to a single node
   - Creates a simplified "shape graph"
   - Similar to the Mapper algorithm in topological data analysis

4. **Analyze the shape graph**
   - Geodesic distances: How far apart are different states?
   - Temporal connectivity matrix (TCM): Recurrence structure
   - Graph topology: Reveals the structure of the attractor


#### Shape Graph
- **Nodes** = Groups of time points with similar states
- **Edges** = Transitions between groups
- For Lotka-Volterra: Should show a **cyclic structure**!

#### Temporal Connectivity Matrix (TCM)
- Shows shortest path length between all time point pairs
- Diagonal patterns = periodic behavior
- For periodic systems: repeated diagonal bands

#### Geodesic Distances
- Shortest paths between nodes on shape graph
- Captures "dynamical distance" not just Euclidean distance
- Opposite phases of cycle have maximum distance

### Sampling interval (SI)
- Should capture dynamics adequately
- ~10-20 samples per oscillation period is good
- Too sparse: miss dynamics; Too dense: slow computation

## Comparison with Neural Data

| Aspect | Lotka-Volterra | Neural (gen_sim.m) |
|--------|----------------|-------------------|
| State space | 2D (prey, predator) | 66D (brain regions) |
| Dynamics | Limit cycle | Attractor landscape |
| Noise | Optional | Always present |
| Ground truth | Known equations | Unknown |
| Interpretation | Population cycles | Brain states |


- Original tmapper repository methods (Mengsen Zhang)
- Lotka, A.J. (1925). Elements of Physical Biology
- Volterra, V. (1926). Variazioni e fluttuazioni del numero d'individui in specie animali conviventi

