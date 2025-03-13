penalty = 30

Properties = {
    'mw': {
        'LLD': 320,
        'LUD': 420,
        'LLA': 170,
        'LUA': 570,
        'N': 5
    },
    'clogp': {
        'LLD': 2,
        'LUD': 3,
        'LLA': 0,
        'LUA': 6,
        'N': 0.1
    },
    'tpsa': {
        'LLD': 40,
        'LUD': 60,
        'LLA': 0,
        'LUA': 100,
        'N': 1
    }
}


def meets_all_requirements(**kwargs):
    for prop, value in kwargs.items():
        if prop not in Properties:
            raise ValueError(f"Propriedade desconhecida: {prop}")
        limits = Properties[prop]
        if not (limits['LLD'] < value < limits['LUD']):
            return False
    return True


def fitness(mw, clogp, tpsa):

    if meets_all_requirements(**locals()):
        return 100

    score = 100

    # Penalty for mw
    if mw < Properties['mw']['LLD']:
      if mw < Properties['mw']['LLA']:
        score -= penalty
      else:        
        score -= (Properties['mw']['LLD'] - mw) / Properties['mw']['N']

    if mw > Properties['mw']['LUD']:
      if mw > Properties['mw']['LUA']:
        score -= penalty
      else:
        score -= (mw - Properties['mw']['LUD']) / Properties['mw']['N']
    
    # Penalty for clogp
    if clogp < Properties['clogp']['LLD']:
      score -= (Properties['clogp']['LLD'] - clogp) / Properties['clogp']['N']

    if clogp > Properties['clogp']['LUD']:
      if clogp >= Properties['clogp']['LUA']:
        score -= penalty
      else:
        score -= (clogp - Properties['clogp']['LUD']) / Properties['clogp']['N']    
    
    # Penalty for tpsa
    if tpsa < Properties['tpsa']['LLD']:
      score -= (Properties['tpsa']['LLD'] - tpsa) / Properties['tpsa']['N']

    if tpsa > Properties['tpsa']['LUD']:
      if tpsa > Properties['tpsa']['LUA']:
        score -= penalty
      else:  
        score -= (tpsa - Properties['tpsa']['LUD']) / Properties['tpsa']['N']
    
    return score
