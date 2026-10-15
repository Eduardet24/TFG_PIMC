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

Despues de cada simulacion, ejecuta `python3 Script_Armonico1D_basico.py` para
generar sus graficas y resumen. Luego mueve los datos brutos y esas salidas a
una carpeta `N_<beads>`. Para una simulacion con `acceptance_rate=50`, usa el
subdirectorio `AR50`. Ejecuta solo el bloque que corresponda a la simulacion;
no ejecutes ambos porque el primer `mv` retira los archivos de origen:

```bash
mkdir -p Resultados/1DArmonico/AR0
mv outrd.dat oute.dat Resultados/1DArmonico/AR0/
mv Resultados/1DArmonico/ejecucion_actual/* Resultados/1DArmonico/AR0/

# O bien, para acceptance_rate=50
mkdir -p Resultados/1DArmonico/AR50/N_10
mv outrd.dat oute.dat Resultados/1DArmonico/AR50/N_10/
mv Resultados/1DArmonico/ejecucion_actual/* Resultados/1DArmonico/AR50/N_10/
```

Una vez archivadas las corridas, ejecuta:

```bash
python3 Script_Armonico1D.py
```

El script analiza los datos brutos de cada carpeta, genera las graficas de
densidad y residuos junto a cada corrida y actualiza las tablas y graficas de
convergencia para la serie de `acceptance_rate=0`, la serie `AR50` y la
comparacion conjunta. No copies las graficas ni el CSV de una corrida a otra:
se regeneran a partir de los datos archivados.
