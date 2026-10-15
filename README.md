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

Con `numpy` y `matplotlib` instalados, ejecutar desde este directorio:

```bash
python3 Script_Armonico1D.py
```

El script guarda los resultados en:

```text
Resultados/1DArmonico/
```

Los archivos de resultados y las graficas se sobrescriben en cada ejecucion.

## Estudio de convergencia con el numero de beads

El script toma el `Npop` de `MC.cfg` y analiza solo los archivos actuales
`outrd.dat` y `oute.dat` del directorio de ejecucion. No lee datos de las
carpetas `N_<beads>`. Despues de ejecutar `./MC` y el script de analisis, puedes
archivar manualmente esa ejecucion antes de iniciar la siguiente. Por ejemplo,
para guardar una ejecucion con 10 beads:

```bash
mkdir -p Resultados/1DArmonico/N_10
cp outrd.dat oute.dat Resultados/1DArmonico/*.png Resultados/1DArmonico/convergencia_beads_1D.csv Resultados/1DArmonico/N_10/
```

La tabla y la grafica de convergencia generadas en cada ejecucion contienen
solo el punto de ese `Npop`, comparado con el valor analitico. Las carpetas
archivadas sirven para conservar manualmente cada ejecucion; el script no las
consulta ni combina sus datos.
