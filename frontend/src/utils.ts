export const getLocalToken = () => localStorage.getItem('token');

export const saveLocalToken = (token: string) => localStorage.setItem('token', token);

export const removeLocalToken = () => localStorage.removeItem('token');

export enum ChannelColor {
  Grayscale = '',
  Red = 'r',
  Green = 'g',
  Blue = 'b',
  Yellow = 'y',
  Cyan = 'c',
  Magenta = 'm'
}

export function convertIndexToColor(index: number): string {
  let color = ChannelColor.Grayscale;
  switch (index) {
    case 0: {
      color = ChannelColor.Red;
      break;
    }
    case 1: {
      color = ChannelColor.Green;
      break;
    }
    case 2: {
      color = ChannelColor.Blue;
      break;
    }
    case 3: {
      color = ChannelColor.Yellow;
      break;
    }
    case 4: {
      color = ChannelColor.Cyan;
      break;
    }
    case 5: {
      color = ChannelColor.Magenta;
      break;
    }
    default: {
      color = ChannelColor.Grayscale;
      break;
    }
  }
  return color;
}

export function convertColorToIndex(color: string): number {
  let index = -1;
  switch (color) {
    case ChannelColor.Red: {
      index = 0;
      break;
    }
    case ChannelColor.Green: {
      index = 1;
      break;
    }
    case ChannelColor.Blue: {
      index = 2;
      break;
    }
    case ChannelColor.Yellow: {
      index = 3;
      break;
    }
    case ChannelColor.Cyan: {
      index = 4;
      break;
    }
    case ChannelColor.Magenta: {
      index = 5;
      break;
    }
    default: {
      index = -1;
      break;
    }
  }
  return index;
}
