from idaes.models.properties.modular_properties.eos.ceos import Cubic



class NRTL(Cubic):
    @staticmethod
    def common(b, pobj):




        # Add activity coefficient parameters as necessary
        if self.config.valid_phase == ("Liq", "Vap") or self.config.valid_phase == (
            "Vap",
            "Liq",
        ):
            # NRTL Model specific variables (values to be fixed by user
            # or need to be estimated based on VLE data)
            # See documentation for suggested or typical values.
            self.alpha = Var(
                self.component_list,
                self.component_list,
                initialize=0.3,
                doc="Non-randomness parameter for NRTL model",
            )

            self.tau = Var(
                self.component_list,
                self.component_list,
                initialize=1.0,
                doc="Binary interaction parameter " "for NRTL model",
                )
        else:
            raise Exception("This property package only supports liquid and vapor phases")



        pass
    @staticmethod
    def build_parameters(b):
        pass

    @staticmethod
    def build_critical_properties(b):
        pass



    @staticmethod
    def _make_NRTL_eq(blk):
        # This code is originally from the IDAES  activity_coeff_prop_pack.py
        # NRTL model variables
        blk.Gij_coeff = Var(
            blk.params.component_list,
            blk.params.component_list,
            initialize=1.0,
            doc="Gij coefficient for use in NRTL model ",
        )

        blk.activity_coeff_comp = Var(
            blk.params.component_list,
            initialize=1.0,
            doc="Activity coefficient of component",
        )

        blk.A = Var(
            blk.params.component_list,
            initialize=1.0,
            doc="Intermediate variable to compute activity" " coefficient",
        )

        blk.B = Var(
            blk.params.component_list,
            initialize=1.0,
            doc="Intermediate variable to compute activity" " coefficient",
        )

        def rule_Gij_coeff(self, i, j):
            # i,j component
            if i != j:
                return self.Gij_coeff[i, j] == exp(
                    -self.params.alpha[i, j] * self.params.tau[i, j]
                )
            else:
                self.Gij_coeff[i, j].fix(1)
                return Constraint.Skip

        blk.eq_Gij_coeff = Constraint(
            blk.params.component_list, blk.params.component_list, rule=rule_Gij_coeff
        )

        # First sum part in the NRTL equation
        def rule_A(self, i):
            value_1 = sum(
                self.mole_frac_phase_comp["Liq", j]
                * self.params.tau[j, i]
                * self.Gij_coeff[j, i]
                for j in self.params.component_list
            )
            value_2 = sum(
                self.mole_frac_phase_comp["Liq", k] * self.Gij_coeff[k, i]
                for k in self.params.component_list
            )
            return self.A[i] == value_1 / value_2

        blk.eq_A = Constraint(blk.params.component_list, rule=rule_A)

        # Second sum part in the NRTL equation
        def rule_B(self, i):
            value = sum(
                (
                    self.mole_frac_phase_comp["Liq", j]
                    * self.Gij_coeff[i, j]
                    / sum(
                        self.mole_frac_phase_comp["Liq", k] * self.Gij_coeff[k, j]
                        for k in self.params.component_list
                    )
                )
                * (
                    self.params.tau[i, j]
                    - sum(
                        self.mole_frac_phase_comp["Liq", m]
                        * self.params.tau[m, j]
                        * self.Gij_coeff[m, j]
                        for m in self.params.component_list
                    )
                    / sum(
                        self.mole_frac_phase_comp["Liq", k] * self.Gij_coeff[k, j]
                        for k in self.params.component_list
                    )
                )
                for j in self.params.component_list
            )
            return self.B[i] == value

        blk.eq_B = Constraint(blk.params.component_list, rule=rule_B)

        # Activity coefficient using NRTL
        def rule_activity_coeff(self, i):
            return log(self.activity_coeff_comp[i]) == self.A[i] + self.B[i]

        blk.eq_activity_coeff = Constraint(
            blk.params.component_list, rule=rule_activity_coeff
        )
    
    @staticmethod
    def _make_flash_eq(blk):
        # This code is originally from the IDAES  activity_coeff_prop_pack.py
        # and modified to suit the generic property package
        if blk.params.config.state_vars == "FTPz":
            # Total mole balance
            def rule_total_mass_balance(self):
                return (
                    self.flow_mol_phase["Liq"] + self.flow_mol_phase["Vap"]
                    == self.flow_mol
                )

            blk.eq_total = Constraint(rule=rule_total_mass_balance)

            # Component mole balance
            def rule_comp_mass_balance(self, i):
                return self.flow_mol * self.mole_frac_comp[i] == (
                    self.flow_mol_phase["Liq"] * self.mole_frac_phase_comp["Liq", i]
                    + self.flow_mol_phase["Vap"] * self.mole_frac_phase_comp["Vap", i]
                )

            blk.eq_comp = Constraint(
                blk.params.component_list, rule=rule_comp_mass_balance
            )

            # sum of mole fractions constraint (sum(x_i)-sum(y_i)=0)
            def rule_mole_frac(self):
                return (
                    sum(
                        self.mole_frac_phase_comp["Liq", i]
                        for i in self.params.component_list
                    )
                    - sum(
                        self.mole_frac_phase_comp["Vap", i]
                        for i in self.params.component_list
                    )
                    == 0
                )

            blk.eq_sum_mol_frac = Constraint(rule=rule_mole_frac)

            if blk.config.defined_state is False:
                # applied at outlet only as complete state information unknown
                blk.eq_mol_frac_out = Constraint(
                    expr=sum(blk.mole_frac_comp[i] for i in blk.params.component_list)
                    == 1
                )
        else:

            def rule_comp_mass_balance(self, i):
                return (
                    self.flow_mol_comp[i]
                    == self.flow_mol_phase_comp["Liq", i]
                    + self.flow_mol_phase_comp["Vap", i]
                )

            blk.eq_comp = Constraint(
                blk.params.component_list, rule=rule_comp_mass_balance
            )

            def rule_mole_frac(self, p, i):
                return (
                    self.mole_frac_phase_comp[p, i]
                    * sum(
                        self.flow_mol_phase_comp[p, i]
                        for i in self.params.component_list
                    )
                    == self.flow_mol_phase_comp[p, i]
                )

            blk.eq_mole_frac = Constraint(
                blk.params._phase_component_set, rule=rule_mole_frac
            )

        # Smooth Flash Formulation

        # Please refer to Burgard et al., "A Smooth, Square Flash
        # Formulation for Equation Oriented Flowsheet Optimization",
        # Computer Aided Chemical Engineering 44, 871-876, 2018.

        blk._temperature_equilibrium = Var(
            initialize=blk.temperature.value,
            doc="Temperature for calculating " "phase equilibrium",
            units=pyunits.K,
        )

        blk._t1 = Var(
            initialize=blk.temperature.value,
            doc="Intermediate temperature for calculating "
            "the equilibrium temperature",
            units=pyunits.K,
        )

        blk.eps_1 = Param(
            default=0.01,
            mutable=True,
            doc="Smoothing parameter for equilibrium " "temperature",
            units=pyunits.K,
        )
        blk.eps_2 = Param(
            default=0.0005,
            mutable=True,
            doc="Smoothing parameter for equilibrium " "temperature",
            units=pyunits.K,
        )

        # Equation #13 in reference cited above
        # Approximation for max(temperature, temperature_bubble)
        def rule_t1(b):
            return b._t1 == 0.5 * (
                b.temperature
                + b.temperature_bubble
                + sqrt((b.temperature - b.temperature_bubble) ** 2 + b.eps_1**2)
            )

        blk._t1_constraint = Constraint(rule=rule_t1)

        # Equation #14 in reference cited above
        # Approximation for min(_t1, temperature_dew)
        # TODO : Add option for supercritical extension
        def rule_teq(b):
            return b._temperature_equilibrium == 0.5 * (
                b._t1
                + b.temperature_dew
                - sqrt((b._t1 - b.temperature_dew) ** 2 + b.eps_2**2)
            )

        blk._teq_constraint = Constraint(rule=rule_teq)

        def rule_phase_eq(self, i):
            return self.fug_phase_comp["Vap", i] == self.fug_phase_comp["Liq", i]

        blk.eq_phase_equilibrium = Constraint(
            blk.params.component_list, rule=rule_phase_eq
        )