#!/bin/sh

rm -f ../tmp/*

echo "GPS Position Calculation Program: start `date`"

../bin/GPSCALPOS >> ../tmp/list

echo "GPS Position Calulation Program: end `date`"
