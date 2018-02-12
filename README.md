# Pathways

We perform a large-scale analysis of *disease pathways* in the human *interactome* to better understand connectivity and higher-order network structure of disease pathways.

Discovering disease pathways, which can be defined as sets of proteins associated with a given disease, is an important problem that has the potential to provide clinically actionable insights for disease diagnosis, prognosis, and treatment.

![Examples of disease pathways](/images/disease-pathway-examples.png)

Here we study the PPI network structure of disease pathways. We find that 90% of pathways do not correspond to single well-connected components in the PPI network. Instead, proteins associated with a single disease tend to form many separate connected components/regions in the network. We then evaluate state-of-the-art disease pathway discovery methods and show that their performance is especially poor on diseases with disconnected pathways.

![Disease pathway discovery](/images/disease-pathway-discovery.png)

We conclude that network connectivity structure alone may not be sufficient for disease pathway discovery. However, we show that higher-order network structures, such as small subgraphs of the pathway, provide a promising direction for the development of new methods.

Please check the [project website](http://snap.stanford.edu/pathways/) for more details. 
  
## Citing

If you find disease pathway analysis useful for your research, please consider citing:

    @inproceedings{agrawal2018,
      author={Agrawal, Monica and Zitnik, Marinka and Leskovec, Jure},
      title = {Large-scale Analysis of Disease Pathways in the Human Interactome},
      year = {2018},
      booktitle = {Pacific Symposium on Biocomputing},
      volume = {23},
      pages = {111-122}
    }

## Miscellaneous

Please send any questions you might have about the code and/or the 
algorithm to <marinka@cs.stanford.edu>.

## License

Decagon is licensed under the MIT License.
