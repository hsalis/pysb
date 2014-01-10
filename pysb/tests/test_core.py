from pysb.testing import *
from pysb.core import *
from functools import partial
from copy import deepcopy
from pysb.core import MonomerPattern

def test_component_names_valid():
    for name in 'a', 'B', 'AbC', 'dEf', '_', '_7', '__a01b__999x_x___':
        c = Component(name, _export=False)
        eq_(c.name, name)

def test_component_names_invalid():
    for name in 'a!', '!B', 'A!bC~`\\', '_!', '_!7', '__a01b  999x_x___!':
        assert_raises(InvalidComponentNameError, Component, name, _export=False)

def test_monomer():
    sites = ['x', 'y', 'z']
    states = {'y': ['foo', 'bar', 'baz'], 'x': ['e']}
    A = Monomer('A', sites, states, _export=False)
    assert_equal(A.sites, sites)
    assert_equal(A.site_states, states)
    assert_equal(type(A()), MonomerPattern)

    assert_raises(ValueError, Monomer, 'A', 'x', _export=False)
    assert_raises(Exception, Monomer, 'A', 'x', 'x', _export=False)
    assert_raises(Exception, Monomer, 'A', ['x'], {'y': ['a']}, _export=False)
    assert_raises(Exception, Monomer, 'A', ['x'], {'x': [1]}, _export=False)

    # Test equality operator
    assert_equal_symmetric(A, A) # Should equal itself
    A_copy = deepcopy(A)
    assert_equal_symmetric(A, A_copy) # Should equal its copy
    A_copy.sites=['y', 'x', 'z']
    assert_equal_symmetric(A, A_copy) # Site order doesn't matter
    A_copy2 = deepcopy(A)
    A_copy2.site_states['y'] = ['foo']
    assert_not_equal_symmetric(A, A_copy2) # Lists contain different members
    A_copy2.site_states['y'] = ['foo', 'baz', 'bar']
    assert_equal_symmetric(A, A_copy2) # State order shouldn't matter
    B = Monomer('B', ['x','y'], _export=False)
    assert_not_equal_symmetric(A, B) # B lacks site 'z', so not equal

def test_unhashable():
    """Mutable objects with __eq__ methods should not be hashable."""
    A = Monomer('A', [], _export=False)
    assert_raises(TypeError, hash, A) # Monomer
    assert_raises(TypeError, hash, A()) # MonomerPattern
    P = Parameter('P', 0, _export=False)
    assert_raises(TypeError, hash, P) # Parameter
    C = Compartment('C', dimension=3, _export=False)
    assert_raises(TypeError, hash, C) # Compartment

def test_monomer_pattern():
    sites = ['x', 'y', 'z']
    states = {'y': ['foo', 'bar', 'baz'], 'x': ['e']}
    m = Monomer('A', sites, states, _export=False)
    mp = MonomerPattern(m, {'x':ANY}, Compartment("Joe", _export=False))
    o  = deepcopy(mp)
    assert_equal(mp, o)
    o.compartment = Compartment("Jim", _export=False)
    ok_(mp != o)
    o  = deepcopy(mp)
    o.site_conditions = {'x':'e'}
    ok_(mp != o)
    o  = deepcopy(mp)
    mp.monomer = Monomer('B', sites, states, _export=False)
    ok_(mp != o)

def test_compartment():
    c = Compartment("Joe", _export=False)
    o = deepcopy(c)
    ok_(c == o)
    o.name = "Jim"
    ok_(c != o)
    o = deepcopy(c)
    o.size = 2
    ok_(c != o)

def test_parameter():
    p = Parameter("a", 2.3, _export=False)
    o = deepcopy(p)
    ok_(o==p)
    o.value = 2.4
    ok_(o!=p)
    o = deepcopy(p)
    o.name = "b"
    ok_(o!=p)

@with_model
def test_monomer_model():
    Monomer('A', ['x','y'])
    ok_(A in model.monomers)
    ok_(A in model.all_components())
    ok_(A not in model.all_components() - model.monomers)


@with_model
def test_initial():
    Monomer('A', ['s'])
    Parameter('A_0')
    Initial(A(s=None), A_0)
    assert_raises_iice = partial(assert_raises, InvalidInitialConditionError,
                                 Initial)
    assert_raises_iice('not a complexpattern', A_0)
    assert_raises_iice(A(), A_0)
    assert_raises_iice(A(s=None), A_0)
    assert_raises_iice(MatchOnce(A(s=None)), A_0)

@with_model
def test_model_pickle():
    import pickle
    A = Monomer('A', _export=False)
    B = Monomer('B', ['x', 'y'], {'x': ['e', 'f']}, _export=False)
    k = Parameter('k', 1.0, _export=False)
    r = Rule('r', A() + B(x='e', y=WILD) >> A() % B(x='f', y=None), k,
             _export=False)
    o = Observable('o', A() % B(), _export=False)
    e = Expression('e', k * o, _export=False)
    c = Compartment('c', None, 3, k, _export=False)
    for comp in [A, B, k, r, o, e, c]:
        model.add_component(comp)
    model.add_component(c)
    Initial(A() ** c, k)
    assert_equal(len(model.all_components()), 7)
    model2 = pickle.loads(pickle.dumps(model))
    check_model_against_component_list(model, model2.all_components())

@with_model
def test_compartment_initial_error():
    Monomer('A', ['s'])
    Parameter('A_0', 2.0)
    c1 = Compartment("C1")
    c2 = Compartment("C2")
    Initial(A(s=None)**c1, A_0)
    Initial(A(s=None)**c2, A_0)
    
