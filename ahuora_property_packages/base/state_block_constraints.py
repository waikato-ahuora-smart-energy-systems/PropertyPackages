from pyomo.environ import Block, Constraint
from pyomo.core.base.expression import ScalarExpression, Expression, _GeneralExpressionData, ExpressionData
from pyomo.core.base.var import ScalarVar, _GeneralVarData, VarData, IndexedVar, Var


class StateBlockConstraints:

    """
    Extended state blocks should inherit from this class.
    Adds additional expressions and constraints to a state block.
    """

    def build(blk, *args):
        blk.constraints = Block()
        blk.add_extra_expressions()


    def constrain(blk, name: str, value: float) -> Constraint | Var | None:
        """constrain a component by name to a value"""
        var = getattr(blk, name)
        return blk.constrain_component(var, value)


    def constrain_component(blk, component: Var | Expression, value: float) -> Constraint | Var | None:
        """
        Constrain a component to a value
        """
        if isinstance(component, ScalarExpression):
            c = Constraint(expr=component == value)
            c.defining_state_var = True
            blk.constraints.add_component(component.local_name, c)
            return c
        elif type(component) in (ScalarVar, _GeneralVarData, VarData, IndexedVar):
            component.fix(value)
            return component
        elif type(component) in (_GeneralExpressionData, ExpressionData):
            # allowed, but we don't need to fix it (eg. mole_frac_comp in helmholtz)
            return None
        else:
            raise Exception(
                f"Component {component} is not a Var or Expression: {type(component)}"
            )


    def add_extra_expressions(blk):
        """
        IDAES state blocks don't support all the properties
        we need, so we add some extra expressions here.
        
        This method can be overridden in a subclass to add
        additional expressions specific to the property package.
        """
        if not hasattr(blk, "enth_mass"):
            blk.add_component("enth_mass", Expression(expr=(blk.flow_mol * blk.enth_mol) / blk.flow_mass))
        if not hasattr(blk, "entr_mass"):
            blk.add_component("entr_mass", Expression(expr=(blk.flow_mol * blk.entr_mol) / blk.flow_mass))
        if not hasattr(blk, "entr_mol"):
            blk.add_component("entr_mol", Expression(expr=(blk.flow_mol * blk.entr_mass) / blk.flow_mass))
        if not hasattr(blk, "total_energy_flow"):
            blk.add_component("total_energy_flow", Expression(expr=blk.flow_mass * blk.enth_mass))
