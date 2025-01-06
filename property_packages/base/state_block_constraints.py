from pyomo.environ import Block, Constraint
from pyomo.core.base.expression import ScalarExpression, Expression, _GeneralExpressionData, ExpressionData
from pyomo.core.base.var import ScalarVar, _GeneralVarData, VarData, IndexedVar, Var

from property_packages.utils.add_extra_expressions import add_extra_expressions


class StateBlockConstraints:

    """
    Extended state blocks should inherit from this class.
    Adds additional expressions and constraints to a state block.
    """

    def build(blk, *args):
        add_extra_expressions(blk)
        blk.constraints = Block()


    def constrain(blk, name: str, value: float) -> Constraint | Var | None:
        """constrain a component by name to a value"""
        var = getattr(blk, name)
        return blk.constrain_component(var, value)


    def constrain_component(blk, component: Var | Expression, value: float) -> Constraint | Var | None:
        """
        Constrain a component to a value
        """
        if type(component) == ScalarExpression:
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
