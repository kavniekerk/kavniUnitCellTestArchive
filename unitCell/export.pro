Group{
Domain_Dirichlet += Region[750000]; 
Domain_Dirichlet += Region[760000]; 
Domain_Dirichlet += Region[770000]; 
Domain_Dirichlet += Region[780000]; 
Domain_Dirichlet += Region[780052]; 
Domain_Dirichlet += Region[780053]; 
DomainCC_Ele += Region[780054]; 
DomainCC_Ele += Region[780055]; 
}
Include "MaterialDatabase.pro";
Function {
epsr[Region[780054]] = Air_epsilonr;
epsr[Region[780055]] = Copper_epsilonr;
}
Constraint { { Name ElectricScalarPotential; Case { { Region Region[750000]; Value -400; } { Region Region[760000]; Value -400; } { Region Region[770000]; Value -400; } { Region Region[780000]; Value -400; } { Region Region[780052]; Value 0; } { Region Region[780053]; Value -2000; } } } }
Constraint { { Name GlobalElectricPotential; Case { } } }
Constraint { { Name GlobalElectricCharge; Case { } } }
Include "Electrostatics.pro";
