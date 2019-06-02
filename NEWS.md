Metrosamp changelog
===================

Version 0.3
-----------

* Preserve attributes of `p0` when generating proposals.  Log-posterior
  functions that rely on attributes (e.g. `names`) can therefore
  rely on those attributes being present.
  
* The `names` attribute, if any, is also copied to `colnames` for the
  samples matrix.
