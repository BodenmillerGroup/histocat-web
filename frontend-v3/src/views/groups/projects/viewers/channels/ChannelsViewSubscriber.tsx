import React from "react";
import TitleInfo from "vitessce/components/TitleInfo";
import { ChannelsView } from "./ChannelsView";
import { CoordinationScopes } from "vitessce/types";

type ChannelsViewSubscriberProps = {
  coordinationScopes: CoordinationScopes;
  removeGridComponent(event: any): void;
};

export function ChannelsViewSubscriber(props: ChannelsViewSubscriberProps) {
  const { removeGridComponent } = props;
  return (
    <TitleInfo
      title="Channels"
      removeGridComponent={removeGridComponent}
      isScroll={true}
      theme="dark"
      isReady={true}
    >
      <ChannelsView />
    </TitleInfo>
  );
}
