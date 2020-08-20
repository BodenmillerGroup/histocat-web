<template>
  <v-card tile class="ma-6 pa-1">
    <v-card-title>
      <v-row no-gutters>
        <v-col>
          <h5 class="text-h5">{{ group.name }}</h5>
          <span v-if="group.url" class="text-caption"><v-icon small>mdi-link</v-icon> {{ group.url }}</span>
        </v-col>
      </v-row>
    </v-card-title>
    <v-card-text v-if="group.description">
      {{ group.description }}
    </v-card-text>
    <v-card-text v-if="group.tags && group.tags.length > 0">
      <v-chip :key="item" v-for="item in group.tags" label small class="mr-1">
        <v-icon small left>mdi-tag-outline</v-icon>
        {{ item }}
      </v-chip>
    </v-card-text>
    <v-card-actions>
      <v-btn
        v-if="isMember || user.is_admin"
        color="primary"
        :to="{ name: 'main-group', params: { groupId: group.id } }"
      >
        Open
      </v-btn>
      <v-btn v-if="!isMember && group.is_open" color="primary" @click="joinGroup()">
        Join
      </v-btn>
      <v-spacer />
      <v-tooltip bottom>
        <template v-slot:activator="{ on }">
          <v-btn icon v-on="on" :to="{ name: 'main-groups-edit', params: { groupId: group.id } }">
            <v-icon>mdi-pencil</v-icon>
          </v-btn>
        </template>
        <span>Edit group</span>
      </v-tooltip>
      <v-tooltip bottom>
        <template v-slot:activator="{ on }">
          <v-btn icon color="secondary" v-on="on" @click="deleteGroup(group.id)">
            <v-icon>mdi-delete-outline</v-icon>
          </v-btn>
        </template>
        <span>Delete group</span>
      </v-tooltip>
    </v-card-actions>
  </v-card>
</template>

<script lang="ts">
import { Component, Prop, Vue } from "vue-property-decorator";
import { groupModule } from "@/modules/group";
import { IUserProfile } from "@/modules/user/models";
import { IGroup } from "@/modules/group/models";

@Component
export default class GroupCard extends Vue {
  readonly groupContext = groupModule.context(this.$store);

  @Prop({
    type: Object,
    required: true,
  })
  readonly user!: IUserProfile;

  @Prop({
    type: Object,
    required: true,
  })
  readonly group!: IGroup;

  get userIds() {
    return this.group.members.filter((item) => item.is_active).map((member) => member.user_id);
  }

  get isMember() {
    return this.userIds.includes(this.user.id);
  }

  async joinGroup() {
    await this.groupContext.actions.joinGroup(this.group.id);
  }

  async deleteGroup() {
    if (self.confirm("Are you sure you want to delete the group?")) {
      await this.groupContext.actions.deleteGroup(this.group.id);
    }
  }
}
</script>
