import pytest
from FCIDUMP_tools.permute_FCIDUMP import expand_perm_multielec, convert_perm_rep, convert_fcidump_idx

class TestExpandPermMultielec:
    @pytest.mark.parametrize("perm,multiplier,expected", [
        # From original test_basic_expansion
        ((1, 4, 3, 2), 2, (1, 2, 7, 8, 5, 6, 3, 4)),
        # From original test_no_expansion  
        ((4, 6, 1, 3, 2, 5), 1, (4, 6, 1, 3, 2, 5)),
        # From original test_single_element_expansion
        ((1,), 5, (1, 2, 3, 4, 5)),
        # Additional test cases
        ((1, 2), 2, (1, 2, 3, 4)),
        ((2, 1), 2, (3, 4, 1, 2)),
        ((1, 2, 3), 1, (1, 2, 3)),
        ((1,), 3, (1, 2, 3)),
    ])
    def test_expand_perm_multielec(self, perm, multiplier, expected):
        result = expand_perm_multielec(perm, multiplier)
        assert result == expected


class TestConvertPermRep:
    @pytest.mark.parametrize("perm,expected", [
        ((1, 2, 3, 4, 5, 6, 7, 8), (1, 2, 3, 4, 5, 6, 7, 8)),
        ((1, 3, 2, 5, 4, 7, 6, 8), (1, 3, 2, 5, 4, 7, 6, 8)),
        ((5, 1, 6, 2, 7, 8, 3, 4), (2, 4, 7, 8, 1, 3, 5, 6))
    ])
    def test_conversion(self, perm, expected):
        result = convert_perm_rep(perm)
        assert result == expected

# TODO: Introduce a test for convert_fcidump_idx