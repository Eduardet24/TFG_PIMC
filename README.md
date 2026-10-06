# Simulacion PIMC del oscilador armonico 1D

## Compilacion rapida

Desde este directorio:

```bash
gcc -O2 -fopenmp *.c -lm -o MC
```

Tambien se puede usar el Makefile de Intel, si `icc` esta disponible:

```bash
make -f Makefile.txt
```

Para eliminar los objetos compilados:

```bash
make -f Makefile.txt clean
```

## Ejecucion

El programa lee los parametros de `MC.cfg` y genera los archivos de salida de la simulacion:

```bash
./MC
```

Entre ellos se encuentran `outrd.dat`, con la distribucion radial, y `oute.dat`, con la energia.

## Analisis y graficas

Con `numpy` y `matplotlib` instalados, ejecutar:

```bash
python3 script_python_graficas.py
```

El script guarda los resultados en:

```text
Resultados/1DArmonico/
```

Los archivos de resultados y las graficas se sobrescriben en cada ejecucion.
