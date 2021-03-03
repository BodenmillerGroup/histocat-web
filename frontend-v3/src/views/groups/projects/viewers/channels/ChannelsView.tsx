import styles from "./ChannelsView.module.scss";
import shallow from "zustand/shallow";
import React, { useState } from "react";
import { Button, InputGroup, Intent } from "@blueprintjs/core";
import { isEqual, throttle } from "lodash-es";
import { EditChannelDialog } from "./EditChannelDialog";
import { useProjectsStore } from "modules/projects";
import { SuggestionView } from "components/SuggestionView";
import { ChannelView } from "./ChannelView";
import { IChannel } from "modules/projects/models";

export function ChannelsView() {
  const { activeAcquisition, selectedMetals, setSelectedMetals, getChannelStackImage } = useProjectsStore(
    (state) => ({
      activeAcquisition: state.getActiveAcquisition(),
      selectedMetals: state.selectedMetals,
      setSelectedMetals: state.setSelectedMetals,
      getChannelStackImage: state.getChannelStackImage,
    }),
    shallow
  );
  const [filterValue, setFilterValue] = useState<string>("");
  const [editDialogOpen, setEditDialogOpen] = useState<boolean>(false);
  const [activeItem, setActiveItem] = useState<IChannel | null>(null);

  if (!activeAcquisition) {
    return <SuggestionView title="" content="Please select acquisition" intent={Intent.NONE} />;
  }

  const channels = Object.values(activeAcquisition.channels).sort((a, b) => a.mass - b.mass);
  const filter = filterValue.toLowerCase();
  const data = channels.filter((item) => {
    return item.name.toLowerCase().includes(filter) || item.label.toLowerCase().includes(filter);
  });

  const handleFilterChange = (event: React.FormEvent<HTMLElement>) =>
    setFilterValue((event.target as HTMLInputElement).value);

  const editAction = (item: IChannel) => {
    setActiveItem(item);
    setEditDialogOpen(true);
  };

  const selectAction = async (channel: IChannel, selected: boolean) => {
    const newSelectedMetals = selected
      ? selectedMetals.concat(channel.name)
      : selectedMetals.filter((metal) => metal !== channel.name);
    if (!isEqual(newSelectedMetals, selectedMetals)) {
      setSelectedMetals(newSelectedMetals);
      await getChannelStackImage();
    }
  };

  return (
    <div className={styles.container}>
      <span className={styles.toolbar}>
        <InputGroup
          asyncControl={true}
          leftIcon="filter"
          rightElement={
            filterValue ? <Button icon="cross" minimal={true} onClick={() => setFilterValue("")} /> : undefined
          }
          onChange={throttle(handleFilterChange, 200, { leading: false })}
          placeholder="Filter channels..."
          value={filterValue}
          fill={true}
        />
      </span>
      <div className={styles.scrollable}>
        {data.map((channel) => (
          <ChannelView
            key={channel.name}
            channel={channel}
            selected={selectedMetals.includes(channel.name)}
            onSelect={selectAction}
            onEdit={editAction}
          />
        ))}
      </div>
      {activeItem && (
        <EditChannelDialog
          channel={activeItem}
          isOpen={editDialogOpen}
          handleClose={() => {
            setEditDialogOpen(false);
            setActiveItem(null);
          }}
        />
      )}
    </div>
  );
}
