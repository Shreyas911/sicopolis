.. _mb_account:

Mass balance accounting
***********************

Lorem ipsum...

Parameter ``MB_ACCOUNT`` ...

..
  #define MB_ACCOUNT 1
                        Mass balance accounting by "hidden ablation scheme"
                        (by R. Calov, A. Robinson):
                         0 : Glaciation of all inner points of the
                             domain allowed
                             (prevents accurate accounting of calving
                              near the margin)
                         1 : Outermost inner points of the domain
                             (i=1,IMAX-1, j=1,JMAX-1) not allowed to glaciate
                             (required for accurate accounting of calving
                              near the margin)
