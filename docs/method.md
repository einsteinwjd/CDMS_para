2. Methods
2.1 Statistical Moment Analysis of Aerosol Size Distributions
This study employed statistical moment analysis to characterize the temporal evolution of aerosol size distributions. Particle number size distribution data were collected across multiple size bins ranging from approximately 1 nm to 1000 nm, providing comprehensive coverage of the nucleation, Aitken, and accumulation modes. The statistical moments approach offers a robust mathematical framework for extracting essential distribution parameters while being computationally efficient.

2.1.1 Moment Calculation
For each time point in the measurement series, we calculated statistical moments from zeroth to sixth order. The k-th order moment ($M_k$) of a particle size distribution is defined as:

$$M_k = \sum_{i=1}^{n} n_i D_i^k$$

where $n_i$ represents the number concentration of particles in the i-th size bin with geometric mean diameter $D_i$, and $n$ is the total number of size bins. These moments provide a mathematical foundation for characterizing the complete distribution through a limited set of parameters.

2.1.2 Derived Distribution Parameters
From these statistical moments, we derived several physically meaningful parameters to characterize the aerosol population:

Total number concentration ($N_{total}$): Directly represented by the zeroth moment $$N_{total} = M_0$$

Mean diameter ($D_{mean}$): Calculated as the ratio of the first moment to the zeroth moment $$D_{mean} = \frac{M_1}{M_0}$$

Surface area ($S$): Proportional to the second moment $$S \propto \pi M_2$$

Volume ($V$): Proportional to the third moment $$V \propto \frac{\pi}{6} M_3$$

Effective diameter ($D_{eff}$): Defined as the ratio of the third moment to the second moment $$D_{eff} = \frac{M_3}{M_2}$$

Geometric standard deviation ($\sigma_g$): Calculated using the relationship between moments $$\sigma_g = \sqrt{\frac{M_2 \cdot M_0}{M_1^2}}$$

These parameters provide comprehensive insights into the size distribution characteristics, including central tendency, dispersion, and higher-order features.

2.2 Distribution Reconstruction and Verification
To verify the validity of the moment-based approach, we reconstructed particle size distributions using the calculated moments and compared them with the original measurement data. Assuming a lognormal distribution, which has been widely used to describe atmospheric aerosol populations, the reconstructed number size distribution function is expressed as:

$$\frac{dN}{d\log D} = \frac{N_{total}}{\sqrt{2\pi}\log\sigma_g}\exp\left[-\frac{(\log D - \log D_g)^2}{2(\log\sigma_g)^2}\right]$$

where $D_g$ is the geometric mean diameter calculated from:

$$D_g = \exp\left(\log D_{mean} - \frac{1}{2}\log\sigma_g^2\right)$$

The consistency between reconstructed distributions and original data was evaluated at selected time points throughout the measurement period.

2.3 Multimodal Distribution Analysis
Atmospheric aerosol size distributions often exhibit multimodal characteristics due to various formation and transformation processes. To identify and characterize these multimodal structures, we employed higher-order moment analysis combined with specialized fitting techniques.

2.3.1 Modal Structure Identification
We utilized the relationship between higher-order moments to detect multimodality in the size distributions. Specifically, we calculated:

Kurtosis: A measure of the "tailedness" of the distribution $$\text{Kurtosis} = \frac{M_5 \cdot M_1}{M_3^2}$$

Skewness: A measure of the asymmetry of the distribution $$\text{Skewness} = \frac{M_4 \cdot M_1}{M_3 \cdot M_2}$$

Distributions with kurtosis values exceeding a threshold (determined empirically to be approximately 4.2 for our dataset) were flagged as potential candidates for multimodal analysis.

2.3.2 Multi-Peak Distribution Fitting
For distributions exhibiting multimodal characteristics, we implemented an adaptive multi-peak lognormal distribution fitting algorithm. The algorithm employs a simulated annealing-like approach with multiple initial parameter sets to optimize the fitting process. For a distribution with N peaks, the fitting model is expressed as:

$$\frac{dN}{d\log D} = \sum_{i=1}^{N} a_i \exp\left[-\frac{(\log D - b_i)^2}{2c_i^2}\right]$$

where $a_i$, $b_i$, and $c_i$ represent the amplitude, logarithmic center position, and width parameter of the i-th peak, respectively.

Our algorithm includes specialized techniques for handling overlapping peaks, including:

Peak number estimation using multi-scale smoothing analysis
Enhanced parameter initialization for overlapping regions
Quality assessment that considers both mathematical fit and physical plausibility

2.4 Temporal Evolution Analysis
To analyze the temporal evolution of aerosol size distributions, we constructed time series of the calculated moments and derived parameters. Additionally, we developed heatmap visualizations to represent the complete size distribution dynamics over time.

For the heatmap analysis, we arranged the size distribution data in a two-dimensional matrix where rows represent time points and columns represent particle size bins. This visualization method allows for intuitive identification of nucleation events, growth processes, and modal shifts in the aerosol population.

The time evolution of modal parameters extracted from the multi-peak fitting was analyzed to identify patterns in aerosol dynamics, including:

Formation and growth rates of nucleation mode particles
Transfer between modes due to coagulation and condensation
Removal processes affecting different size ranges
This comprehensive analytical framework provides a robust methodology for characterizing complex aerosol size distribution data, supporting detailed investigations of atmospheric aerosol dynamics and transformation processes