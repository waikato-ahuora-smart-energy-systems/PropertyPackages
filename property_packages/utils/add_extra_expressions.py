from pyomo.environ import Block, Expression, value, units

def add_extra_expressions(sb: Block):
    """
    idaes state blocks don't support all the properties we need, so we add some extra expressions here
    """
    if not hasattr(sb, "enth_mass"):
        sb.add_component("enth_mass", Expression(expr=(sb.flow_mol * sb.enth_mol) / sb.flow_mass))
    if not hasattr(sb, "entr_mass"):
        sb.add_component("entr_mass", Expression(expr=(sb.flow_mol * sb.entr_mol) / sb.flow_mass))
    if not hasattr(sb, "entr_mol"):
        sb.add_component("entr_mol", Expression(expr=(sb.flow_mol * sb.entr_mass) / sb.flow_mass))
    if not hasattr(sb, "total_energy_flow"):
        sb.add_component("total_energy_flow", Expression(expr=sb.flow_mass * sb.enth_mass))