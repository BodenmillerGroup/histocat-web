import styles from "./ChannelView.module.scss";
import { IChannel } from "modules/projects/models";
import { Button, Checkbox } from "@blueprintjs/core";

type ChannelViewProps = {
  channel: IChannel;
  selected: boolean;
  onSelect: (channel: IChannel, selected: boolean) => void;
  onEdit: (channel: IChannel) => void;
};

export function ChannelView(props: ChannelViewProps) {
  return (
    <span className={styles.container}>
      <Checkbox
        checked={props.selected}
        onChange={(event) => props.onSelect(props.channel, (event.target as HTMLInputElement).checked)}
      />
      <span>{props.channel.name}</span>
      <span className="bp3-text-small bp3-text-muted">{props.channel.customLabel}</span>
      <Button minimal={true} small={true} icon="edit" onClick={() => props.onEdit(props.channel)} />
    </span>
  );
}
