from amp.descriptor.gaussian import Gaussian
from amp.utilities import hash_images

def get_fingerprints(images):
    descriptor=Gaussian()
    images=hash_images(images, ordered=True)
    descriptor.calculate_fingerprints(images)
