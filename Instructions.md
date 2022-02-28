# statisticalmodelingofturbulence


1. Run autospec_BendatPiersol.m in Matlab to see how autospectrum is extracted from a wind database. 
2. A sample wind data file (20141031-ModInt-B02.mat) is provided, and it will be processed by the Matlab code by assigning mnum = 9. To run this code, wppwerrornoA.m, filenames.mat, psdperturbseries.m and PSDofRAW.m are required.
3. The CNN model is MultiparaCNN-KFold.ipynb. To run this code, you will also need 5paracnn5crawdisjointPSDlogspacefreq.m. 
4. Similarly, the LSTM model is LSTM-Kfold-minicurve.ipynb. To run this code, you will also need logmldata500.mat.
5. The LSTM model can be compared with the perturbation series results using the Matlab file compare_nnmodel_minicurve.m. Download the mincurve folder to run this code.
