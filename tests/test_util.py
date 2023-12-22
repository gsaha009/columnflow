# coding: utf-8


__all__ = ["UtilTest"]

import unittest

from columnflow.util import (
    create_random_name, maybe_import, MockModule, DotDict, Derivable,
    safe_div, test_float, test_int, is_regex, is_pattern, pattern_matcher,
)


class UtilTest(unittest.TestCase):

    def test_create_random_name(self):
        self.assertIsInstance(create_random_name(), str)

    def test_maybe_import(self):
        # non-existing package
        mock_module = maybe_import("not_existing_package")
        self.assertIsInstance(mock_module, MockModule)
        with self.assertRaises(Exception):
            mock_module()

        with self.assertRaises(ImportError):
            _ = maybe_import("not_existing_package", force=True)

        # existing package
        self.assertEqual(unittest, maybe_import("unittest"))

    def test_save_div(self):
        self.assertEqual(safe_div(10, 2), 5)
        self.assertEqual(safe_div(10, 0), 0)

    def test_test_int_test_float(self):
        for test_number in (test_int, test_float):
            self.assertIsInstance(test_number(1), bool)
            self.assertIsInstance(test_number("some_string"), bool)
            self.assertTrue(test_number(1))
            self.assertTrue(test_number(1.5))
            self.assertTrue(test_number(-1.5))
            self.assertFalse(test_number("some_string"))
            self.assertFalse(test_number(1j))
            self.assertFalse(test_number([1, 2]))

    def test_is_regex(self):
        self.assertTrue(is_regex(r"^foo\d+.*$"))
        self.assertFalse(is_regex(r"^no$atEnd"))
        self.assertFalse(is_regex(r"no^atStart$"))

    def test_is_pattern(self):
        self.assertTrue(is_pattern("foo*"))
        self.assertTrue(is_pattern("bar?"))
        self.assertFalse(is_pattern("not_a_pattern"))

    def test_pattern_matcher(self):
        matcher = pattern_matcher("foo*")
        self.assertTrue(matcher("foo"))
        self.assertTrue(matcher("foo123"))
        self.assertFalse(matcher("xfoo"))
        self.assertFalse(matcher("bar123"))

        matcher = pattern_matcher(r"^foo\d+.*$")
        self.assertTrue(matcher("foo1"))
        self.assertTrue(matcher("foo12x3"))
        self.assertFalse(matcher("foo"))
        self.assertFalse(matcher("foox1"))

        matcher = pattern_matcher(("foo*", "*bar"), mode=any)
        self.assertTrue(matcher("foo123"))
        self.assertTrue(matcher("123bar"))

        matcher = pattern_matcher(("foo*", "*bar"), mode=all)
        self.assertFalse(matcher("foo123"))
        self.assertFalse(matcher("123bar"))
        self.assertTrue(matcher("foo123bar"))

    def test_DotDict(self):
        my_dict = DotDict({
            "A": 1,
            "B": 2,
        })
        self.assertEqual(my_dict.A, 1)
        self.assertEqual(my_dict.get("B", 3), 2)
        self.assertEqual(my_dict.get("C", 3), 3)
        with self.assertRaises(AttributeError):
            my_dict.C

        my_dict = DotDict.wrap({
            "A": 1,
            "B": 2,
            "C": {"D": 3},
        })
        self.assertEqual(my_dict.A, 1)
        self.assertEqual(my_dict.C.D, 3)
        my_dict = DotDict.wrap({
            "A": 1,
            "B": 2,
            "C": [3, 4],
        })
        self.assertEqual(my_dict.C, [3, 4])

    def test_Derivable(self):
        class MyClass(Derivable):
            A = 1
            B = 2

        derived_class = MyClass.derive("derived", cls_dict={
            "A": 3,
            "C": 4,
            "dummy_method": lambda N: 3 * N,
        })
        deep_derived_class = derived_class.derive("deep_derived")

        self.assertEqual(derived_class.cls_name, "derived")
        self.assertEqual(derived_class.A, 3)
        self.assertEqual(derived_class.B, 2)
        self.assertEqual(derived_class.C, 4)
        with self.assertRaises(AttributeError):
            derived_class.D
        with self.assertRaises(AttributeError):
            MyClass.C

        self.assertEqual(MyClass._subclasses, {"derived": derived_class})
        self.assertTrue(MyClass.has_cls("derived"))
        self.assertFalse(MyClass.has_cls("deep_derived", deep=False))
        self.assertTrue(MyClass.has_cls("deep_derived", deep=True))
        self.assertEqual(MyClass.get_cls("deep_derived", deep=True), deep_derived_class)
        self.assertFalse(MyClass.has_cls("not_derived"))
