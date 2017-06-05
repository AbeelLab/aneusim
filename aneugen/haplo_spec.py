import ast
import operator as op

from aneugen.haplogen import MutationType, SPEC_MUTTYPE_KEY, DISTANCE_MODELS

# supported operators for math evaluation
operators = {ast.Add: op.add, ast.Sub: op.sub, ast.Mult: op.mul,
             ast.Div: op.truediv, ast.Pow: op.pow, ast.BitXor: op.xor,
             ast.USub: op.neg}


def eval_math(expr):
    """
    Evaluate simple math expressions.

    From
    https://stackoverflow.com/questions/2371436/evaluating-a-mathematical-expression-in-a-string

    >>> eval_expr('2^6')
    4
    >>> eval_expr('2**6')
    64
    >>> eval_expr('1 + 2*3**(4^5) / (6 + -7)')
    -5.0
    """
    return _math_eval(ast.parse(expr, mode='eval').body)


def _math_eval(node):
    if isinstance(node, ast.Num):  # <number>
        return node.n
    elif isinstance(node, ast.BinOp):  # <left> <operator> <right>
        return operators[type(node.op)](_math_eval(node.left),
                                        _math_eval(node.right))
    elif isinstance(node, ast.UnaryOp):  # <operator> <operand> e.g., -1
        return operators[type(node.op)](_math_eval(node.operand))
    else:
        raise TypeError(node)


def get_distance_model(chromosome_spec, mut_type: MutationType):
    mut_type_key = SPEC_MUTTYPE_KEY[mut_type]
    modelname_key = "{}.model".format(mut_type_key)

    # First check if this mutation type is disabled by checking for boolean
    # value
    try:
        modelname = chromosome_spec.getboolean(modelname_key)
        if not modelname:
            return None
    except ValueError:
        pass

    # Otherwise try to get the actual distance model
    modelname = chromosome_spec.get(modelname_key)

    if modelname not in DISTANCE_MODELS:
        raise KeyError("Invalid distance model '{}'".format(modelname))

    cls = DISTANCE_MODELS[modelname]
    return cls.from_spec(chromosome_spec, mut_type)


def get_dosage(chromosome_spec):
    if 'ploidy' not in chromosome_spec:
        raise KeyError("No ploidy information found")

    ploidy = chromosome_spec.getint('ploidy')

    dosage_key = 'dosage_dist'
    if dosage_key not in chromosome_spec:
        alt_key = 'dosage_dist.{}'.format(ploidy)

        if alt_key not in chromosome_spec:
            raise KeyError("No dosage information found.")

        dosage_key = alt_key

    dosages = chromosome_spec.get(dosage_key).split(',')

    # Evaluate simple mathematical expressions, so they can enter fractions as
    # dosage values
    dosages = list(map(eval_math, (v.strip() for v in dosages)))

    if len(dosages) != ploidy:
        raise ValueError("Not enough dosage values. Got {} values, expected "
                         "{} values.".format(len(dosages), ploidy))

    return dosages
