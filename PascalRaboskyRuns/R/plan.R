plan <- drake_plan(
  phy = ape::rcoal(30),
  hisse_out = DoSingleRun(phy)
)
