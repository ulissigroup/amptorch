import random
from amptorch.active_learning.bootstrap import bootstrap_ensemble


def test_bootstrap():
    random.seed(0)

    data = [1]
    ensembles = 3
    ensemble_sets, parent_dataset = bootstrap_ensemble(
        parent_dataset=data, n_ensembles=ensembles
    )

    assert len(ensemble_sets) == ensembles, "Incorrect number of ensembles returned"
    for ensemble in ensemble_sets:
        assert (
            ensemble == data
        ), "Ensemble sets with parent_dataset length 1 is incorrect!"
    assert data == parent_dataset, "Parent dataset is incorrect for initial length 1"

    new_data = 3
    # manually computed bootstrap ensemble with same random seed
    new_parent = [1, 3]
    e_1 = [1, 3]
    e_2 = [1, 3]
    e_3 = [1, 3]

    ensemble_sets, parent_dataset = bootstrap_ensemble(
        parent_dataset=data,
        resampled_set=ensemble_sets,
        new_data=new_data,
        n_ensembles=ensembles,
    )

    assert len(ensemble_sets) == ensembles, "Incorrect number of ensembles returned"
    assert ensemble_sets[0] == e_1, "Ensemble set is incorrect!"
    assert ensemble_sets[1] == e_2, "Ensemble set is incorrect!"
    assert ensemble_sets[2] == e_3, "Ensemble set is incorrect!"

    assert (
        new_parent == parent_dataset
    ), "Parent dataset is incorrect when adding new data!"

    data = [1, 2, 5, 6, 9]
    ensembles = 3
    ensemble_sets, parent_dataset = bootstrap_ensemble(
        parent_dataset=data, n_ensembles=ensembles
    )

    #TODO write test for initial parent_datset > 1
