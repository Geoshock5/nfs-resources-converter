import { ChangeDetectionStrategy, Component } from '@angular/core';
import { GuiComponentInterface } from '../../gui-component.interface';

@Component({
  selector: 'app-palette.block-ui',
  templateUrl: './palette.block-ui.component.html',
  styleUrls: ['./palette.block-ui.component.scss'],
  changeDetection: ChangeDetectionStrategy.OnPush,
})
export class PaletteBlockUiComponent implements GuiComponentInterface {

  resourceData: ReadData | null = null;
  name: string = '';

  constructor() { }

  lpad(str: string, padString: string, length: number) {
    while (str.length < length)
      str = padString + str;
    return str;
  }
}