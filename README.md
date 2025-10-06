# TCRAD - Tripleclouds Radiance Model

This document last updated 6 October 2025

Robin Hogan <r.j.hogan@ecmwf.int>

## INTRODUCTION

This package contains the C++ version of the 1D Triplecluds Radiance
Model (TCRAD) that calculates thermal-infrared and microwave radiances
in cloudy scenes, suitable for use in data assimilation. The
"Tripleclouds" method for representing cloud structure divides each
cloudy layer into three regions, one clear and two cloudy with
different optical depth. The exchange of radiation between the regions
of adjacent layers is governed by exponential-random overlap rules.

The software makes use of the
[Adept](https://github.com/rjhogan/Adept-2) array and automatic
differentiation library, through which the tangent-linear and adjoint
interfaces to the radiance calculations are provided.

A Fortran interface is provided, using Fortran-2003 interoperability
features to allow Fortran assumed-shape arrays to be passed in.


## PACKAGE OVERVIEW

The subdirectories are as follows:

- `src` - the TCRAD souce code

- `include` - Fortran include file


## TO COMPILE

First create the configure script:

  autoreconf -i

At ECMWF type

  module swap gcc/14.2.0
  ./configure_ecmwf.sh
  make

You can edit the script to change the configuration.

Outside ECMWF you first need to install the Adept library from
https://github.com/rjhogan/Adept-2 and call the configure script with
"--with-adept=/your/adept/install/directory".


## LICENCE

(C) Copyright 2014- ECMWF.

This software is licensed under the terms of the Apache Licence Version 2.0
which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

In applying this licence, ECMWF does not waive the privileges and immunities
granted to it by virtue of its status as an intergovernmental organisation
nor does it submit to any jurisdiction.
Copyright statements are given in the file NOTICE.



## CONTACT

Please email Robin Hogan <r.j.hogan@ecmwf.int> with any queries or bug
fixes, but note that ECMWF does not commit to providing support to
users of this software.
