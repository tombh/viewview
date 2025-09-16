# pip2 install pillow

import logging
import sys
import os
import shutil
from optparse import OptionParser
from PIL import Image

logger = logging.getLogger(__name__)

def is_folder_a_tiletree(folder):
    is_tiletree = True
    for f in os.listdir(folder):

        if(os.path.isdir(os.path.join(folder,f))):
            try:
                _i = int(f)
                if(_i > 18): # Google maps only have zoom levels to 19
                    is_tiletree = False
                    break
            except ValueError, e:
                is_tiletree = False
                break
    return is_tiletree

ACCEPTABLE_IMAGE_EXTS = set(['.png', '.gif'])

def remove_alpha(img):
    return img
    img = img.convert("RGBA")
    datas = img.getdata()

    newData = []
    for item in datas:
        if item[3] > 0 and item[3] < 255:
            newData.append((item[0], item[1], item[2], 0))
        else:
            newData.append(item)

    img.putdata(newData)
    return img

def merge_tiletrees(src_folder, dst_root_folder,
                    dry_run=False, action='copy'):
    if not os.path.exists(dst_root_folder):
        os.mkdir(dst_root_folder)

    if not is_folder_a_tiletree(src_folder):
        raise RuntimeError("Source folder doesn't seem to be a tiletree." \
                               " It should contain folders with zoom levels," \
                               " ie 9/ 10/ 11/ ...\n\n")

    if action == 'copy':
        action_f = shutil.copy
    elif action == 'move':
        action_f = shutil.move
    else:
        raise RuntimeError("Invalid action: %s" % action)

    for root, dirs, files in os.walk(src_folder):
        canonical_dir = os.path.relpath(root, src_folder)
        dest_dir = os.path.join(dst_root_folder, canonical_dir)

        if not os.path.exists(dest_dir):
            logger.info("Making dir %s " % dest_dir)
            os.mkdir(dest_dir)

        for f in files:
            name, ext = os.path.splitext(f)
            if ext.lower() not in ACCEPTABLE_IMAGE_EXTS:
                logger.debug("Skipping file '%s'" % f)
                continue
            src_file = os.path.join(root,f)
            canonical_file = os.path.join(canonical_dir, f)

            dest_file = os.path.join(dest_dir,f)
            if os.path.exists(dest_file):
                if dry_run:
                    logger.debug("Dry run: Would have merged file %s -> %s"
                                 % (src_file, dest_file))
                else:
                    logger.debug("Merging file %s -> %s"
                                 % (src_file, dest_file))
                    img1 = Image.open(src_file)
                    img2 = Image.open(dest_file)
                    img1 = Image.composite(img1, remove_alpha(img2), remove_alpha(img1))
                    img1.save(dest_file)
            else:
                if dry_run:
                    logger.debug("Dry run: Would have [%s] from %s -> %s"
                                 % (action, src_file, dest_file))
                else:
                    logger.debug("%s from %s -> %s"
                                 % (action, src_file, dest_file))
                    action_f(src_file,dest_dir)

if __name__ == "__main__":
    parser = OptionParser()
    parser.add_option("-d", action="store_true", dest="dry_run", default=False,
                      help="Don't perform any copy/moves")
    parser.add_option("-v", action="store_true", dest="verbose", default=False,
                      help="Verbose output")
    parser.add_option("--move", action="store_true", dest="move", default=False,
                      help="Move images from the source folder instead of" \
                          " copying. Can be much faster, but use with care!")
    (options, args) = parser.parse_args()

    log_level = logging.DEBUG if options.verbose else logging.INFO
    logging.basicConfig(format='%(message)s', level=log_level)

    action = 'move' if options.move else 'copy'
    src_folder = os.path.normpath(args[0])
    dst_root_folder = os.path.normpath(args[1])

    merge_tiletrees(src_folder, dst_root_folder,
                    dry_run=options.dry_run, action=action)
