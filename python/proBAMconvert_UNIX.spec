# -*- mode: python -*-

block_cipher = None


a = Analysis(['proBAM.py'],
             pathex=['/home/vladie/PycharmProjects/proBAMconvert/python'],
             binaries=[],
             datas=[('proBAMconvert_logo.gif', '.')],
             hiddenimports=['pysam.libctabixproxies'],
             hookspath=[],
             runtime_hooks=[],
             excludes=['scipy', 'matplotlib', 'IPython', 'sqlalchemy', 'PyQt4'],
             win_no_prefer_redirects=False,
             win_private_assemblies=False,
             cipher=block_cipher)
pyz = PYZ(a.pure, a.zipped_data,
             cipher=block_cipher)
exe = EXE(pyz,
          a.scripts,
          a.binaries,
          a.zipfiles,
          a.datas,
          name='proBAMconvert_UNIX',
          debug=False,
          strip=False,
          upx=True,
          console=False )
