import os
import importlib
import unittest

class TestImportFunctions(unittest.TestCase):

    def test_imports(self):
        package_dir = 'morph_utils'
        files = os.listdir(package_dir)
        for file in files:
            print(file)
            if file.endswith('.py') and file != '__init__.py':
                module_name = os.path.splitext(file)[0]
                full_module_path = f'{package_dir}.{module_name}'

                with self.subTest(module=full_module_path):
                    try:
                        module = importlib.import_module(full_module_path)
                        self.assertIsNotNone(module, f'Failed to import module: {full_module_path}')

                        # Check if functions can be accessed
                        for attribute_name in dir(module):
                            attribute = getattr(module, attribute_name)
                            if callable(attribute):
                                self.assertIsNotNone(attribute, f'Failed to import function: {attribute_name} in module: {full_module_path}')

                    except Exception as e:
                        self.fail(f'Error importing module {full_module_path}: {str(e)}')

if __name__ == '__main__':
    unittest.main()
