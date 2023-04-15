import pytest

from lib import resource_calculator as rc


@pytest.fixture
def queue_set():
    queues = {"small": ["q1", "q2", "q3"], "large": ["q4", "q5"], "huge": ["q6"]}
    return queues


@pytest.mark.parametrize(
    "selected_queue,expected_queues",
    [("small", ["q1", "q2", "q3"]), ("large", ["q4", "q5"]), ("huge", ["q6"])],
)
def test_select_queue(queue_set, selected_queue, expected_queues):
    """
    Test select_queue when the user selection is valid
    """
    res = [rc.select_queue(selected_queue, queue_set) for i in range(100)]
    assert set(res) == set(expected_queues)


def test_select_queue_invalid_selection(queue_set):
    """
    Test select_queue when the user selection is invalid
    """
    with pytest.raises(
        ValueError, match=r"Configured queue set does not match anything in user resource config.*"
    ):
        rc.select_queue("fake_queue", queue_set)
