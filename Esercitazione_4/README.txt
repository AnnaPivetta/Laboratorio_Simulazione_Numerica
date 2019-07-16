Il Jupyter notebook è contenuto nella cartella MolecularDynamics_NVE.
La cartella contiene il codice e alcune cartelle dove sono salvati i file necessari per fare i grafici (sul Jupyter notebook sono caricati i file da queste cartelle, perchè l'esecuzione del codice richiede del tempo).
Per seguire il codice:
mettersi nella cartella MolecularDynamics_NVE
mettere la variabile restart=0 nel file input.dat (in questa fase si può fare anche un solo blocco da 1000 passi)
compilare (make)
eseguire (./MolDyn_NVE)
mettere poi restart=1 nel file di input.dat 
eseguire nuovamente fino a che non si è soddisfatti della temperatura (il cui valore istantaneo è stampato a video. Possono andare bene tre o quattro esecuzioni)
aumentare il numero di blocchi e il numero di passi da file di input ed eseguire nuovamente il codice
