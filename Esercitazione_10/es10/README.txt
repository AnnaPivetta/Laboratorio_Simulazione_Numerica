Il Jupyter-notebook contiene grafici realizzati con dati salvati in alcune cartelle.
Per generare i dati:
andare nella cartella es10.1
selezionare i valori desiderati sul file di input. Per la variabile forma scegliere zero (cerchio) o uno (quadrato). T è la temperatura iniziale.
compilare (make)
eseguire (./es10)
Per il parallelo:
andare della cartella parallelo
selezionare i valori desiderati nel file di input
compilare (make)
eseguire(mpiexec -np NumeroDiCore ./es10)
