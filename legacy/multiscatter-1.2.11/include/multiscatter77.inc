! multiscatter77.inc - Fortran 77 include file for multiple scattering algorithms
!
!  Copyright (C) 2010-2011  University of Reading
!
!  Author: Robin Hogan <r.j.hogan@reading.ac.uk>
!
!  Licensed under the Apache License, Version 2.0 (the "License");
!  you may not use this file except in compliance with the License.
!  You may obtain a copy of the License at
!
!      http://www.apache.org/licenses/LICENSE-2.0
!
!  Unless required by applicable law or agreed to in writing, software
!  distributed under the License is distributed on an "AS IS" BASIS,
!  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!  See the License for the specific language governing permissions and
!  limitations under the License.
!

! Possible status codes returned the functions below
      integer MS_SUCCESS                 ! Function executed successfully
      integer MS_FAILURE                 ! A generic error occurred
      integer MS_MEMORY_ALLOCATION_ERROR ! Error occurred allocating memory
      integer MS_INPUT_ERROR             ! Problem with the input or config data
      integer MS_INVALID_CONTEXT         ! Specified context does not exist
      parameter(MS_SUCCESS = 0)
      parameter(MS_FAILURE = 1)
      parameter(MS_MEMORY_ALLOCATION_ERROR = 2)
      parameter(MS_INPUT_ERROR = 3)
      parameter(MS_INVALID_CONTEXT = 4)

! The ms_set_algorithms subroutine (see below) takes three arguments,
! the context, a code indicating the small-angle algorithm and a code
! indicating the wide-angle algorithm.  Codes for small-angle
! algorithm (PVC = photon variance-covariance method):
      integer MS_SINGLE_AND_SMALL_ANGLE_NONE
      integer MS_SINGLE_SCATTERING          
      integer MS_SMALL_ANGLE_PVC_ORIGINAL    ! Hogan (Applied Optics, 2006)
      integer MS_SMALL_ANGLE_PVC_FAST        ! Hogan (J Atmos Sci, 2008)
      integer MS_SMALL_ANGLE_PVC_EXPLICIT    ! Eloranta-like but PVC framework

      parameter(MS_SINGLE_AND_SMALL_ANGLE_NONE= 0)
      parameter(MS_SINGLE_SCATTERING          = 1)
      parameter(MS_SMALL_ANGLE_PVC_ORIGINAL   = 2)
      parameter(MS_SMALL_ANGLE_PVC_FAST       = 3)
      parameter(MS_SMALL_ANGLE_PVC_EXPLICIT   = 4)


!   Codes for wide-angle algorithm (TDTS = time-dependent two-stream method):
      integer MS_WIDE_ANGLE_NONE
      integer MS_WIDE_ANGLE_TDTS_FORWARD_LOBE   ! Hogan and Battaglia (J Atmos Sci, 2008)
      integer MS_WIDE_ANGLE_TDTS_NO_FORWARD_LOBE! Radar version of TDTS

      parameter(MS_WIDE_ANGLE_NONE                  = 0)
      parameter(MS_WIDE_ANGLE_TDTS_FORWARD_LOBE     = 1)
      parameter(MS_WIDE_ANGLE_TDTS_NO_FORWARD_LOBE  = 2)


! FUNCTIONS TO CONFIGURE THE MULTIPLE SCATTERING CALCULATION

! The configuration information for a multiple scattering calculations
! is stored in a "context", and any number of contexts may be used by
! the same program.  Typically a program would obtain a new context ID
! with the function "ms_new_context", configure it with a set of
! additional subroutine calls, then call one of the multiscatter*
! functions to run the algorithm on a profile of input data.  The
! multiscatter* functions may be called again with the same context,
! with or without changing the configuration information (such as the
! altitude of the instrument which in the case of a satellite instrument
! may change for each profile) in between.

! To see how to call these, look at the Fortran 2003 include file
! "multiscatter.inc"

      integer ms_new_context
      external ms_new_context

      integer multiscatter
      external multiscatter

      integer multiscatter_hsrl
      external multiscatter_hsrl

      integer multiscatter_small_angle_only
      external multiscatter_small_angle_only

      integer multiscatter_small_angle_only_hsrl
      external multiscatter_small_angle_only_hsrl

      integer multiscatter_simple
      external multiscatter_simple

      integer multiscatter_AD
      external multiscatter_AD

      integer multiscatter_hsrl_AD
      external multiscatter_hsrl_AD

      integer multiscatter_small_angle_only_AD
      external multiscatter_small_angle_only_AD

      integer multiscatter_small_angle_only_hsrl_AD
      external multiscatter_small_angle_only_hsrl_AD

      integer multiscatter_simple_AD
      external multiscatter_simple_AD
