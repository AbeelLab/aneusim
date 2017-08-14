from abc import ABCMeta, abstractmethod

import numpy
from scipy.stats import lognorm, expon


class BaseModel(metaclass=ABCMeta):
    @classmethod
    @abstractmethod
    def from_spec(cls, chromosome_spec, prefix):
        """Return a new instance based on information in a "spec" file (read
        using a `configparser.ConfigParser`, and available the variables
        available in `chromosome_spec`."""

    @abstractmethod
    def __call__(self) -> int:
        """Return a new random value according to the model."""


class LogNormalModel(BaseModel):
    def __init__(self, mean: float, sigma: float, loc: float=0.0):
        self.mean = mean
        self.sigma = sigma
        self.loc = loc

        self.rv = lognorm(s=sigma, scale=numpy.exp(mean), loc=loc)

    @classmethod
    def from_spec(cls, chromosome_spec, prefix):
        mean_key = "{}.mean".format(prefix)
        sigma_key = "{}.sigma".format(prefix)
        loc_key = "{}.loc".format(prefix)

        if mean_key not in chromosome_spec or sigma_key not in chromosome_spec:
            raise KeyError("Chromosome specification does not contain the "
                           "parameters for the log-normal distribution."
                           " Keys checked: {}, {}".format(mean_key, sigma_key))

        return cls(chromosome_spec.getfloat(mean_key),
                   chromosome_spec.getfloat(sigma_key),
                   chromosome_spec.getfloat(loc_key, 0.0))

    def __call__(self) -> int:
        return int(round(self.rv.rvs()))


class ExponentialModel(BaseModel):
    def __init__(self, lambda_: float):
        self.rv = expon(scale=1/lambda_)

    @classmethod
    def from_spec(cls, chromosome_spec, prefix):
        rate_key = "{}.lambda".format(prefix)

        if rate_key not in chromosome_spec:
            raise KeyError("Chromosome specification does not contain the rate"
                           " parameter for the exponential distribution. "
                           "Keys checked: {}".format(rate_key))

        return cls(chromosome_spec.getfloat(rate_key))

    def __call__(self):
        return int(round(self.rv.rvs()))


class FixedValueModel(BaseModel):
    def __init__(self, value: int):
        self.value = value

    @classmethod
    def from_spec(cls, chromosome_spec, prefix):
        value_key = "{}.value".format(prefix)

        if value_key not in chromosome_spec:
            raise KeyError("Chromosome specification does not contain the "
                           "value parameter for the fixed value model."
                           " Keys checked: {}".format(value_key))

        return cls(chromosome_spec.getint(value_key))

    def __call__(self):
        return self.value


MODELS = {
    'lognormal': LogNormalModel,
    'exponential': ExponentialModel,
    'fixed': FixedValueModel
}


def get_model(chromosome_spec, prefix):
    model_key = "{}.model".format(prefix)

    if model_key not in chromosome_spec:
        raise KeyError("Chromosome specification does not provide a model name"
                       " (key '{}' not found).".format(model_key))

    try:
        modelname = chromosome_spec.getboolean(model_key)
        if not modelname:
            return None
    except ValueError:
        pass

    modelname = chromosome_spec[model_key]

    if modelname not in MODELS:
        raise ValueError(
            "'{}' is not a valid model for key '{}'. Available models: "
            "{}.".format(modelname, model_key, ", ".join(MODELS.keys()))
        )

    return MODELS[modelname].from_spec(chromosome_spec, prefix)
