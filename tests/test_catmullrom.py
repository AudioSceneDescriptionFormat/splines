import pytest
from splines import CatmullRom


def test_repeated_vertices_with_alpha():
    CatmullRom([[1, 1], [1, 1]], alpha=0)
    with pytest.raises(ValueError, match='(?i)repeated vertices.*alpha'):
        CatmullRom([[1, 1], [1, 1]], alpha=1)
