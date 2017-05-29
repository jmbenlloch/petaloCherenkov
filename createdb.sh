#!/bin/bash

time python createdb.py 0 4 &
time python createdb.py 1 4 & 
time python createdb.py 2 4 & 
time python createdb.py 3 4 &
